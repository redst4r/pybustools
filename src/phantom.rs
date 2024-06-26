use std::{collections::HashMap, hash::Hash};
use bustools::{consistent_genes::{find_consistent, GeneId, InconsistentResolution, MappingMode, EC}, consistent_transcripts::{find_consistent_transcripts, MappingResultTranscript, TranscriptId}, io::{BusFolder, BusReader, BusRecord, BusWriter, BusWriterPlain}, iterators::{CbUmiGroupIterator, CellGroupIterator}, merger::MultiIterator, utils::int_to_seq};
use bustools_cli::butterfly::{self, CUHistogram};
use pyo3::{prelude::*, types::{PyDict, PyFrozenSet, PyTuple}};
use rand::Rng;
use rustphantompurger::phantompurger::make_fingerprint_histogram;
use std::fmt::Debug;

use crate::get_spinner;
/// Quck howto


#[pyfunction]
/// Creates a histogram of fingerprints for the cell barcodes.
/// 
/// For each sample/busfile, records how many UMIs [reads] of a particular CB were seen -> Fingerprint
/// Instead of saving every fingerprint (there's as many as theres CBs), just keep a histogram.
/// 
/// Note: We can count the abundance of a CB in terms of reads, or UMIs
/// 
/// Returns:
/// - samplenames
/// - Histogram of fingerprints (measured in UMI abundance)
/// - Histogram of fingerprints (measured in read abundance)
pub (crate) fn phantom_fingerprint_cb(py: Python<'_>, busfolders: HashMap<String, String>) -> PyResult<(PyObject, PyObject, PyObject)> {

    let (samplenames, fp_umi, fp_read ) = _phantom_fingerprint_cb(busfolders);
    let pythondict_umi = PyDict::new(py);
    let pythondict_reads = PyDict::new(py);

    for (k,v) in fp_umi {
        // need to convert it to a frozenSet explicitly (regular sets cant be hashed in python)
        let fset = PyTuple::new(py, k.iter());
        pythondict_umi.set_item(fset, v)?;
    }
    for (k,v) in fp_read {
        // need to convert it to a frozenSet explicitly (regular sets cant be hashed in python)
        let fset = PyTuple::new(py, k.iter());
        pythondict_reads.set_item(fset, v)?;
    }
    Ok((samplenames.to_object(py), pythondict_umi.to_object(py), pythondict_reads.to_object(py)))
}

#[pyfunction]
pub (crate)  fn phantom_fingerprint_cbumi(py: Python<'_>, busfolders: HashMap<String, String>) -> PyResult<(PyObject, PyObject, PyObject)> {

    let (samplenames, fp_umi, fp_read ) = _phantom_fingerprint_cbumi(busfolders);
    let pythondict_umi = PyDict::new(py);
    let pythondict_reads = PyDict::new(py);

    for (k,v) in fp_umi {
        // need to convert it to a frozenSet explicitly (regular sets cant be hashed in python)
        let fset = PyTuple::new(py, k.iter());
        pythondict_umi.set_item(fset, v)?;
    }
    for (k,v) in fp_read {
        // need to convert it to a frozenSet explicitly (regular sets cant be hashed in python)
        let fset = PyTuple::new(py, k.iter());
        pythondict_reads.set_item(fset, v)?;
    }
    Ok((samplenames.to_object(py), pythondict_umi.to_object(py), pythondict_reads.to_object(py)))
}


/// Main meat of `phantom_fingerprint_cb()`, see its docs
fn _phantom_fingerprint_cb(busfolders: HashMap<String, String>) -> (Vec<String>, HashMap<Vec<usize>, usize>, HashMap<Vec<usize>, usize>){
    let mut h = HashMap::new();
    for (k,filename) in busfolders{
        h.insert(k, BusReader::new(&filename).groupby_cb());
    }
    general_phantom_fingerprint(h)
}

/// Main meat of `phantom_fingerprint_cbumi()`, see its docs
fn _phantom_fingerprint_cbumi(busfolders: HashMap<String, String>) -> (Vec<String>, HashMap<Vec<usize>, usize>, HashMap<Vec<usize>, usize>){
    let mut h = HashMap::new();
    for (k,filename) in busfolders{
        h.insert(k, BusReader::new(&filename).groupby_cbumi());
    }
    general_phantom_fingerprint(h)
}


/// Creates a histogram of fingerprints for busfiles ( or rather, grouped iterators over [Vec<BusRecord>]).
/// 
/// For each sample/busfile, records how many UMIs (reads) of a particular group (typically CB or CB/UMI) were seen -> Fingerprint
/// Instead of saving every fingerprint (there's as many as theres CBs), just keep a histogram.
/// 
/// Note: We can count the abundance of a group in terms of reads, or UMIs
/// 
/// Returns:
/// - samplenames
/// - Histogram of fingerprints (measured in UMI abundance)
/// - Histogram of fingerprints (measured in read abundance)
fn general_phantom_fingerprint<I,K>(h: HashMap<String, I>) -> (Vec<String>, HashMap<Vec<usize>, usize>, HashMap<Vec<usize>, usize>)
where
    I: Iterator<Item = (K, Vec<BusRecord>)>,
    K: Ord + Eq + Debug + Copy,
    
{
    let samplenames: Vec<_> = h.keys().cloned().collect();

    let multi_iter = MultiIterator::new(h);

    let bar = get_spinner();
    let interval = 100_000;

    let mut umi_hist : HashMap<Vec<usize>, usize> = HashMap::new();
    let mut read_hist : HashMap<Vec<usize>, usize> = HashMap::new();

    for (i, (_cb, record_dict)) in multi_iter.enumerate() {

        let mut  umi_fp = Vec::new();
        let mut read_fp = Vec::new();

        for s in samplenames.iter(){
            let (numi, nreads) = match record_dict.get(s) {
                Some(records) => {
                    let numi = records.len();
                    let nreads: u32 = records.iter().map(|r| r.COUNT).sum() ;
                    (numi, nreads as usize)
                },
                None => (0,0),
            };
            umi_fp.push(numi);
            read_fp.push(nreads);
        }
        let val = umi_hist.entry(umi_fp).or_insert(0);
        *val += 1;

        let val = read_hist.entry(read_fp).or_insert(0);
        *val += 1;

        if i % interval == 0 {
            bar.inc(interval as u64);
        }
    }
    (samplenames, umi_hist, read_hist)
}

#[pyfunction]
/// Filter busfiles for hopped reads. 
/// If we detect a CB/UMI in more than one busfile, it gets filtered 
/// (i.e. written to `busfiles_removed`) and NOT written to `busfiles_filtered`
pub (crate)  fn rustphantom_filter(py: Python<'_>, 
    busfiles_input: HashMap<String, String>, 
    busfiles_filtered: HashMap<String, String>, 
    busfiles_removed: HashMap<String, String>, 
) -> PyResult<()> {


        let mut folders = HashMap::new();
        for (k,filename) in busfiles_input.iter(){
            folders.insert(k, BusReader::new(filename));
        }


        let mut buswriters: HashMap<String,BusWriterPlain> = busfiles_filtered.iter()
        .map(|(sname, fname)| 
            (
                sname.to_owned(), 
                BusWriterPlain::new(
                    fname, 
                    folders.get(sname).unwrap().get_params().clone()
                )
            ) 
        )
        .collect();

    let mut buswriters_removed: HashMap<String,BusWriterPlain> = busfiles_removed.iter()
        .map(|(sname, fname)| 
            (
                sname.to_owned(), 
                BusWriterPlain::new(
                    fname, 
                    folders.get(sname).unwrap().get_params().clone()
                )
            ) 
        )
        .collect();

    let mut iterators = HashMap::new();
    for (k,filename) in busfiles_input{
        iterators.insert(k, BusReader::new(&filename).groupby_cbumi());
    };

    let mult = MultiIterator::new(iterators);
    let bar = get_spinner();
    let interval = 100_000;
    for (i, (_cbumi, record_dict)) in mult.enumerate() {
        if record_dict.len() > 1 {
            // skipping, writing into removed
            for (sname, records) in record_dict {
                let wri = buswriters_removed.get_mut(&sname).expect("must exist");   
                wri.write_records(&records);
            }     
        } else {
            for (sname, records) in record_dict {
                let wri = buswriters.get_mut(&sname).expect("must exist");   
                wri.write_records(&records);
            }
        }
        if i % interval == 0 {
            bar.inc(interval as u64);
        }       
    };
    Ok(())
}


// SLOW due to EC mapping
#[pyfunction]
pub (crate) fn rustphantom(py: Python<'_>, busfolders: HashMap<String, String>, t2gfile: &str) -> PyResult<(PyObject, PyObject)> {
    
    let mut h = HashMap::new();
    for (k,filename) in busfolders{
        h.insert(k, BusFolder::new(&filename));
    }

    let fingerprint = make_fingerprint_histogram(h, t2gfile);
    
    // let dir = tempdir().unwrap();
    // let file_path = dir.path().join("phantom.bus");
    // let tmpfilename = file_path.to_str().unwrap();

    // fingerprint.to_csv(tmpfilename);

    let samplenames = fingerprint.get_samplenames();
    let histogram = fingerprint.get_histogram();
    // histogram needs to be converted, its keys are lists atm, whihc cant be hashed in python


    let pythondict = PyDict::new(py);

    for (k,v) in histogram {
        // need to convert it to a frozenSet explicitly (regular sets cant be hashed in python)
        let fset = PyTuple::new(py, k.iter());
        pythondict.set_item(fset, v)?;
    }
    Ok((samplenames.to_object(py), pythondict.to_object(py)))
}

