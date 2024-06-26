use std::collections::HashMap;
use bustools::{consistent_genes::{find_consistent, GeneId, InconsistentResolution, MappingMode}, consistent_transcripts::{find_consistent_transcripts, MappingResultTranscript, TranscriptId}, io::BusFolder, iterators::{CbUmiGroupIterator, CellGroupIterator}};
use bustools_cli::butterfly::{self, CUHistogram};
use pyo3::prelude::*;

use crate::{counting::cbumi_overlap, get_spinner};
/// Quck howto
/// ```python
/// import rustpybustools
/// rustpybustools.make_ecs('/home/michi/bus_testing/bus_output_shortest/', '/home/michi/bus_testing/transcripts_to_genes.txt', True)
/// 

/// A Python module implemented in Rust.
#[pymodule]
fn pybustools(_py: Python, m: &PyModule) -> PyResult<()> {
    // m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(cbumi_overlap, m)?)?;
    m.add_function(wrap_pyfunction!(make_ecs_cb, m)?)?;
    m.add_function(wrap_pyfunction!(make_ecs_ec, m)?)?;
    m.add_function(wrap_pyfunction!(make_ecs_across_cb, m)?)?;

    // m.add_function(wrap_pyfunction!(crate::phantom::phantom_fingerprint_cb, m)?)?;
    // m.add_function(wrap_pyfunction!(crate::phantom::phantom_fingerprint_cbumi, m)?)?;
    // m.add_function(wrap_pyfunction!(crate::phantom::rustphantom, m)?)?;
    // m.add_function(wrap_pyfunction!(crate::phantom::rustphantom_filter, m)?)?;

    m.add_function(wrap_pyfunction!(crate::counting::cb_overlap, m)?)?;
    m.add_function(wrap_pyfunction!(crate::counting::detect_cb_frequencies, m)?)?;
    m.add_function(wrap_pyfunction!(crate::counting::kmer_counter_cb, m)?)?;
    m.add_function(wrap_pyfunction!(crate::counting::count_reads_per_EC, m)?)?;
    
    #[pyfn(m)]
    #[pyo3(name = "make_ecs")]
    /// Count the frequency of frequency for the CB/UMI in the specified busfolders.
    /// Essentially a wrapper around `butterfly::make_ecs`.
    /// # Returns
    /// - a HashMap/dictionary, where keys are the number of reads and values are the number of CB/UMI with that nreads.
    /// 
    fn make_ecs(busfile: &str, matrix_file: &str, transcript_file: &str, t2g_file: &str, collapse_ec:bool) -> PyResult<HashMap<usize, usize>> {

        let folder = BusFolder::from_files(busfile, matrix_file, transcript_file);
        
        let mode = if collapse_ec {
            // println!("Making mapper");
            let mapper = folder.make_mapper(t2g_file);
            MappingMode::Gene(mapper, InconsistentResolution::IgnoreInconsistent)
        } else {
            MappingMode::EC(InconsistentResolution::IgnoreInconsistent)
        };
        // println!("Done Making mapper");

        let h = butterfly::make_ecs(&folder.get_busfile(), mode);
        // println!("umis: {}", h.get_numis());
        Ok(h.into())  // conversion into a hashmap/dict
    }

    Ok(())
}


/// Per barcode, creates a CUHistogram: In that cell, how often do you see a UMI with n_reads
#[pyfunction]
fn make_ecs_across_cb(
    busfile: &str, matrix_file: &str, transcript_file: &str,
    // collapse_ec:bool,
) -> PyResult<HashMap<u64, HashMap<usize, usize>>> {

    let folder = BusFolder::from_files(busfile, matrix_file, transcript_file);
    // let mut h_read: HashMap<u64, CUHistogram> = HashMap::new();
    let mut h_read: HashMap<u64, HashMap<usize, usize>> = HashMap::new();
    
    for (_cb, records) in folder.get_iterator().groupby_cb() {
        let mut h_read_cb = HashMap::new();
        for r in records {
            let freq = h_read_cb.entry(r.COUNT as usize).or_insert(0);
            *freq += 1;
        }
        // let cu_umi = CUHistogram::new(h_read_cb);
        h_read.insert(_cb, h_read_cb);
    }
    
    Ok(h_read.into()  // conversion into a hashmap/dict
    ) 

}
/// Counts frequencies of barcodes per cell.
/// THIS IS NOT "PER BARCODE"
/// Returns two dictionaries (each key is a barcode):
/// - the frequency of a barcode in terms of UMIs
/// - the frequency of a barcode in terms of reads
#[pyfunction]
fn make_ecs_cb(
    busfile: &str, matrix_file: &str, transcript_file: &str,
    // collapse_ec:bool,
) -> PyResult<(HashMap<usize, usize>, HashMap<usize, usize>)> {

    let folder = BusFolder::from_files(busfile, matrix_file, transcript_file);
    
    // let mode = if collapse_ec {
    //     // println!("Making mapper");
    //     let mapper = folder.make_mapper(t2g_file);
    //     MappingMode::Gene(mapper, InconsistentResolution::IgnoreInconsistent)
    // } else {
    //     MappingMode::EC(InconsistentResolution::IgnoreInconsistent)
    // };
    // println!("Done Making mapper");
    let mut h_umi: HashMap<usize, usize> = HashMap::new();
    let mut h_read: HashMap<usize, usize> = HashMap::new();

    for (_cb, records) in folder.get_iterator().groupby_cb() {
        let nreads: u32 = records.iter().map(|r| r.COUNT).sum();
        let v = h_read.entry(nreads as usize).or_insert(0);
        *v += 1;

        let numi = records.len();
        let v = h_umi.entry(numi).or_insert(0);
        *v += 1;
    }
    let cu_umi = CUHistogram::new(h_umi);
    let cu_reads = CUHistogram::new(h_read);

    Ok(
        (cu_umi.into(),cu_reads.into())  // conversion into a hashmap/dict
    ) 
}

/// Counts frequencies of barcodes.
/// Returns two dictionaries (each key is a barcode):
/// - the frequency of a barcode in terms of UMIs
/// - the frequency of a barcode in terms of reads
#[pyfunction]
fn make_ecs_ec(busfile: &str, matrix_file: &str, transcript_file: &str, t2g_file: &str, collapse_ec:bool, mapping_mode: &str) -> PyResult<(HashMap<String, HashMap<usize, usize>>, HashMap<usize, usize>, HashMap<usize, usize>)> {
    let folder = BusFolder::from_files(busfile, matrix_file, transcript_file);

    let results = if mapping_mode == "gene" {
        make_ecs_ec_genemapping(folder, t2g_file, collapse_ec)
    } else if mapping_mode == "transcript" {
        make_ecs_ec_transcript_mapping(folder, collapse_ec)
        
    } else {
        panic!("unkown mapping mode; needs to be [gene;transcript]")
    };

    Ok(results)
}

fn make_ecs_ec_transcript_mapping(folder: BusFolder, collapse_ec:bool)
 -> (HashMap<String, HashMap<usize, usize>>, HashMap<usize, usize>, HashMap<usize, usize>) {

    println!("Making mapper");
    let ecmapper = folder.make_mapper_transcript();
    println!("Done Making mapper");

    let mut h_read: HashMap<(TranscriptId, usize), usize> = HashMap::new();
    let mut h_read_multimapped: HashMap<usize, usize> = HashMap::new();
    let mut h_read_inconsistent: HashMap<usize, usize> = HashMap::new();

    let all_transcripts = ecmapper.get_transcript_list();

    // iterate over all CB/UMI, and if the records map to a single gene
    // add an entry of #reads per UMI into `h_read`
    println!("Iterating cb/umi");
    let spinner = get_spinner();
    let mut counter = 0;
    for (_cbumi, records) in folder.get_iterator().groupby_cbumi() {
    
        match find_consistent_transcripts(&records, &ecmapper) {
            MappingResultTranscript::SingleTranscript(geneid) => {
                let nreads: u32 = records.iter().map(|r| r.COUNT).sum();
                let v = h_read.entry((geneid, nreads as usize)).or_insert(0);
                *v+=1;
            },
            MappingResultTranscript::Multimapped(_glist) => {
                let nreads: u32 = records.iter().map(|r| r.COUNT).sum();
                let v = h_read_multimapped.entry(nreads as usize).or_insert(0);
                *v+=1;

            },
            MappingResultTranscript::Inconsistent => {
                let nreads: u32 = records.iter().map(|r| r.COUNT).sum();
                let v = h_read_inconsistent.entry(nreads as usize).or_insert(0);
                *v+=1;

            },
        }
        counter += 1;
        if counter % 1_000_000 == 0{
            spinner.inc(1_000_000)
        }
    }

    // aggregate into a transcript -> CU_histogram map
    let mut final_read: HashMap<TranscriptId, CUHistogram> = HashMap::new();
    for ((tid, freq), count) in h_read {

        if final_read.contains_key(&tid) {

            // update existing hashmap
            let h = final_read.get_mut(&tid).unwrap();
            h.update(freq, count);

        } else {
            let mut inner = HashMap::new();
            inner.insert(freq, count);
            let h = CUHistogram::new(inner);
            final_read.insert(tid, h);

        }
        // let h= final_read.entry(tid).or_insert()
    }

    // convert tid to actual trasncript name
    let mut final_with_string: HashMap<String, HashMap<usize, usize>> = HashMap::new();
    for (tid, cu) in final_read {
        let ix = tid.0 as usize;
        final_with_string.insert(all_transcripts[ix].0.clone(), cu.into());
    }

    (final_with_string, h_read_inconsistent, h_read_multimapped)
}

fn make_ecs_ec_genemapping(folder: BusFolder, t2g_file: &str, collapse_ec:bool
) -> (HashMap<String, HashMap<usize, usize>>, HashMap<usize, usize>, HashMap<usize, usize>) {

    println!("Making mapper");
    let ecmapper = folder.make_mapper(&t2g_file);
    println!("Done Making mapper");
    // let mapping_mode =  if collapse_ec{
    //      MappingMode::Gene(ecmapper, InconsistentResolution::IgnoreInconsistent)
    // } else {
    //     MappingMode::EC(InconsistentResolution::IgnoreInconsistent)
    // };

    // let mapping_mode = MappingMode::Gene(ecmapper, InconsistentResolution::IgnoreInconsistent);

    // let mut h_umi: HashMap<(GeneId, usize), usize> = HashMap::new();
    let mut h_read: HashMap<(GeneId, usize), usize> = HashMap::new();
    let mut h_read_multimapped: HashMap<usize, usize> = HashMap::new();
    let mut h_read_inconsistent: HashMap<usize, usize> = HashMap::new();

    let all_genes = ecmapper.get_gene_list();

    // iterate over all CB/UMI, and if the records map to a single gene
    // add an entry of #reads per UMI into `h_read`
    println!("Iterating cb/umi");
    let spinner = get_spinner();
    let mut counter = 0;
    for (_cbumi, records) in folder.get_iterator().groupby_cbumi() {
    
        match find_consistent(&records, &ecmapper) {
            bustools::consistent_genes::MappingResult::SingleGene(geneid) => {
                // increment that gene + single observation
                // let v = h_umi.entry((geneid, 1)).or_insert(0);
                // *v+=1;

                let nreads: u32 = records.iter().map(|r| r.COUNT).sum();
                let v = h_read.entry((geneid, nreads as usize)).or_insert(0);
                *v+=1;
            },
            bustools::consistent_genes::MappingResult::Multimapped(_glist) => {
                let nreads: u32 = records.iter().map(|r| r.COUNT).sum();
                let v = h_read_multimapped.entry(nreads as usize).or_insert(0);
                *v+=1;

            },
            bustools::consistent_genes::MappingResult::Inconsistent => {
                let nreads: u32 = records.iter().map(|r| r.COUNT).sum();
                let v = h_read_inconsistent.entry(nreads as usize).or_insert(0);
                *v+=1;

            },
        }
        counter += 1;
        if counter % 1_000_000 == 0{
            spinner.inc(1_000_000)
        }
    }

    // aggregate into a gene -> CU_histogram map
    let mut final_read: HashMap<GeneId, CUHistogram> = HashMap::new();
    for ((geneid, freq), count) in h_read {

        if final_read.contains_key(&geneid) {

            // update existing hashmap
            let h = final_read.get_mut(&geneid).unwrap();
            h.update(freq, count);

        } else {
            let mut inner = HashMap::new();
            inner.insert(freq, count);
            let h = CUHistogram::new(inner);
            final_read.insert(geneid, h);

        }
        // let h= final_read.entry(geneid).or_insert()
    }
    // convert gene id to actual gene name
    let mut final_with_string: HashMap<String, HashMap<usize, usize>> = HashMap::new();
    for (geneid, cu) in final_read {
        let ix = geneid.0 as usize;
        final_with_string.insert(all_genes[ix].0.clone(), cu.into());
    }

    (final_with_string, h_read_inconsistent, h_read_multimapped)
}
