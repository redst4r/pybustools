use std::collections::HashMap;
use bustools::{io::{BusFolder, BusReader, BusRecord}, iterators::{CbUmiGroupIterator, CellGroupIterator}, merger::MultiIterator, utils::int_to_seq};
use pyo3::{prelude::*, types::{PyDict, PyFrozenSet}};
use rand::Rng;
use std::fmt::Debug;

use crate::get_spinner;

// #[pyfunction]
// /// Inside a single busfile, how simlilar are the given barcodes in terms of their UMIs.
// /// Uuseful to hunt down possible barcode-duplicates (i.e. a single cell represented by multuiple barcodes)
// pub (crate) fn umi_overlap_in_cb_pairs(py: Python<'_>, busfile: &str, cbs: Vec<u64>) -> PyResult<PyObject> {

//     // let mut results: HashMap<(u64, u64), (usize, usize)> = HashMap::new();
//     let pythondict = PyDict::new(py);
//     let breader1 = BusReader::new(busfile);

//     for pair in breader1.groupby_cb().filter(|(cb, _records)| cbs.contains(&cb)).combinations(2) {
//         let (cb1 , records1) = &pair[0];
//         let (cb2 , records2) = &pair[1];
//         let umis1: HashSet<u64> = records1.iter().map(|r| r.UMI).collect();
//         let umis2: HashSet<u64> = records2.iter().map(|r| r.UMI).collect();
        
//         let n_inter = umis1.intersection(&umis2).count();
//         let n_union = umis1.union(&umis2).count();
//         // results.insert((*cb1, *cb2), (n_inter, n_union));
//         pythondict.set_item((*cb1, *cb2), (n_inter, n_union))?
//     }
//     Ok(pythondict.to_object(py))

// }



#[pyfunction]
/// Calculates the overlap of busfiles in terms of their CB/UMI 
/// (not considering their frequencies, as does `phantom_fingerprint_cb`)
pub (crate) fn cbumi_overlap(py: Python<'_>, busfolders: HashMap<String, String>) -> PyResult<PyObject> {

    let overlap = detect_cbumi_overlap(busfolders);
    let pythondict = PyDict::new(py);

    for (k,v) in overlap {
        // need to convert it to a frozenSet explicitly (regular sets cant be hashed in python)
        let fset = PyFrozenSet::new(py, k.iter())?;
        pythondict.set_item(fset, v)?;
    }
    Ok(pythondict.to_object(py))
}

#[pyfunction]
/// Calculates the overlap of busfiles in terms of their CBs. 
/// (not considering their frequencies, as does `phantom_fingerprint_cb`)
pub (crate) fn cb_overlap(py: Python<'_>, busfolders: HashMap<String, String>) -> PyResult<PyObject> {

    let overlap = detect_cb_overlap(busfolders);
    let pythondict = PyDict::new(py);

    for (k,v) in overlap {
        // need to convert it to a frozenSet explicitly (regular sets cant be hashed in python)
        let fset = PyFrozenSet::new(py, k.iter())?;
        pythondict.set_item(fset, v)?;
    }
    Ok(pythondict.to_object(py))
}


/// Given a set of iterators over busfiles (keyed/grouped by CB, or CB+UMI) count how many
/// keys overlap between busfiles
fn general_detect_overlap<I,K>(h: HashMap<String, I>) -> HashMap<Vec<String>, usize> 
where
    I: Iterator<Item = (K, Vec<BusRecord>)>,
    K: Ord + Eq + Debug + Copy,
{
    let multi_iter = MultiIterator::new(h);
    let mut counter: HashMap<Vec<String>, usize> = HashMap::new();

    let bar = get_spinner();
    let interval = 100_000;

    for (i, (_key, record_dict)) in multi_iter.enumerate() {
        let mut the_set: Vec<String> = record_dict.keys().cloned().collect();
        the_set.sort();
        let val = counter.entry(the_set).or_insert(0);
        *val += 1;

        if i % interval == 0 {
            bar.inc(interval as u64);
        }
    }
    counter

}

#[pyfunction]
pub (crate) fn detect_cb_frequencies(py: Python<'_>, busfolders: HashMap<String, String>, min_clonesize: usize) -> PyResult<PyObject> {

    let frequencies = get_biggest_clones(busfolders, min_clonesize );

    // for (k,v) in fp_umi {
    //     // need to convert it to a frozenSet explicitly (regular sets cant be hashed in python)
    //     let fset = PyTuple::new(py, v.iter());
    //     pythondict_umi.set_item(k, fset)?;
    // }
    // for (k,v) in fp_read {
    //     // need to convert it to a frozenSet explicitly (regular sets cant be hashed in python)
    //     let fset = PyTuple::new(py, v.iter());
    //     pythondict_reads.set_item(k, fset)?;
    // }
    Ok(frequencies.to_object(py))
}

fn valmap_ref<K, V, V2, F>(fun: F, the_map: &HashMap<K, V>) -> HashMap<K, V2>
where
    F: Fn(&V) -> V2,
    K: Eq + std::hash::Hash + Clone,
{
    let r: HashMap<K, V2> = the_map.iter().map(|(k, v)| (k.clone(), fun(v))).collect();
    r
}

fn get_biggest_clones(busfolders: HashMap<String, String>, min_clonesize: usize) -> Vec<HashMap<String, usize>>{
    let mut h = HashMap::new();
    for (k,filename) in busfolders{
        h.insert(k, BusReader::new(&filename).groupby_cb());
    }
    let samplenames: Vec<_> = h.keys().cloned().collect();


    let multi_iter = MultiIterator::new(h);

    let bar = get_spinner();
    let interval = 100_000;

    let mut clonesize_per_sample: Vec<HashMap<String, usize>> = Vec::new();

    for (i, (_cb, record_dict)) in multi_iter.enumerate() {

        let mut n_records = valmap_ref(|x| x.len(), &record_dict);
        if *n_records.values().max().unwrap() > min_clonesize {

            // fill all samples with zero if no freq was observed
            for s in samplenames.iter() {
                if !n_records.contains_key(s) {
                    n_records.insert(s.to_string(), 0);
                }
            }
            n_records.insert("Barcode".to_string(), _cb as usize);
            clonesize_per_sample.push(n_records);
        }
     
        if i % interval == 0 {
            bar.inc(interval as u64);
        }        
    };

    clonesize_per_sample
}

fn detect_cbumi_overlap(busfolders: HashMap<String, String>) -> HashMap<Vec<String>, usize> {

    let mut h = HashMap::new();
    for (k,filename) in busfolders{
        h.insert(k, BusReader::new(&filename).groupby_cbumi());
    }
    general_detect_overlap(h)
}

fn detect_cb_overlap(busfolders: HashMap<String, String>) -> HashMap<Vec<String>, usize> {

    let mut h = HashMap::new();
    for (k,filename) in busfolders{
        h.insert(k, BusReader::new(&filename).groupby_cb());
    }
    general_detect_overlap(h)
}






#[pyfunction]
pub (crate) fn kmer_counter_cb(busfile: &str, kmer_size: usize, fraction: f64) -> (HashMap<String, usize> , HashMap<String, usize>) {
    let b = BusReader::new(busfile);

    let mut rng = rand::thread_rng();

    // let mut kmer_counter:HashMap<String, usize> = HashMap::with_capacity(4_usize.pow(kmer_size as u32));
    let mut kmer_counter2 = HashMap::with_capacity(4_usize.pow(kmer_size as u32));
    let mut base_counter: HashMap<String, usize> = HashMap::with_capacity(4);
    let bar = get_spinner();
    let interval = 100_000;
    for (i, (cb, _records)) in b.groupby_cb().enumerate() {

        // to subsample the cells, takes long otherwise
        let r: f64 = rng.gen();
        if r > fraction {
            continue
        }

        let seq = int_to_seq(cb, 20);

        for i in 0..seq.len() - kmer_size {
            let slice = &seq[i..i+kmer_size];
            match kmer_counter2.get_mut(slice){
                Some(val) => {*val+=1},
                None => {kmer_counter2.insert(slice.to_owned(), 0_usize);},
            };
        }

        for i in 0..seq.len() {
            let slice = &seq[i..i+1];
            match base_counter.get_mut(slice){
                Some(val) => {*val+=1},
                None => {base_counter.insert(slice.to_owned(), 0_usize);},
            };
        }

        if i % interval == 0 {
            bar.inc(interval as u64);
        }        
    }
    (kmer_counter2, base_counter)
}

fn seq_to_kmers(s: &str, kmer_length: usize) -> Vec<String>{
    /*
    highly annoying due to the UTF-8 BS...
    also very slow, creates alot of String instances on the way
     */
    // let kmers: Vec<_> = s.as_bytes().windows(kmer_length).collect();
    let t: Vec<_> = s.chars().collect();
    t.windows(kmer_length).map(String::from_iter).collect()
    // kmers
}



#[pyfunction]
/// Filter busfiles for hopped reads. 
/// If we detect a CB/UMI in more than one busfile, it gets filtered 
/// (i.e. written to `busfiles_removed`) and NOT written to `busfiles_filtered`
pub (crate) fn count_reads_per_ec(_py: Python<'_>, 
    busfolder: &str, 
) -> PyResult<HashMap<usize, usize>> {

    let bus = BusFolder::new(busfolder);
    let mut ec_vector = HashMap::new();

    for (_cbumi, records) in bus.get_iterator().groupby_cbumi() {
        for r in records {
            let v = ec_vector.entry(r.EC as usize).or_insert(0);
            *v += r.COUNT as usize;
        }
    };

    Ok(ec_vector)
}

