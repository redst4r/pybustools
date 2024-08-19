use std::collections::{BTreeSet, HashMap, HashSet};
use bustools::{consistent_genes::{find_consistent, GeneId, Genename, InconsistentResolution, MappingMode, MappingResult, EC}, consistent_transcripts::{find_consistent_transcripts, MappingResultTranscript, TranscriptId}, io::{BusFolder, BusReader}, iterators::{CbUmiGroupIterator, CellGroupIterator}};
use bustools_cli::butterfly::{self, CUHistogram};
use itertools::Itertools;
use pyo3::{prelude::*, types::{PyDict, PyTuple}};
use polars::prelude::*;
use crate::get_spinner;

/// A wrapper of bustools-cli::CUHistogram
/// We need taht to implement a python conversion trait!
pub (crate) struct MyCUHistogram(CUHistogram);

// impl From<CUHistogram> for MyCUHistogram {
//     fn from(value: CUHistogram) -> Self {
//         MyCUHistogram(value)
//     }
// }

/// Convenience, allows us to return CUHistograms directly
impl pyo3::IntoPy<pyo3::PyObject> for MyCUHistogram {
    fn into_py(self, py: pyo3::Python<'_>) -> pyo3::PyObject {
        self.0.get_histogram().to_object( py )
    }
}


#[pyfunction]
/// Count the frequency of frequency for the CB/UMI in the specified busfolders.
/// Essentially a wrapper around `butterfly::make_ecs`.
/// # Returns
/// - a HashMap/dictionary, where keys are the number of reads and values are the number of CB/UMI with that nreads.
/// 
pub (crate) fn make_ecs(busfile: &str, matrix_file: &str, transcript_file: &str, t2g_file: &str, collapse_ec:bool) -> PyResult<MyCUHistogram> {

    let folder = BusFolder::from_files(busfile, matrix_file, transcript_file);
    let mode = if collapse_ec {
        // println!("Making mapper");
        let mapper = folder.make_mapper(t2g_file);
        MappingMode::Gene(mapper, InconsistentResolution::IgnoreInconsistent)
    } else {
        MappingMode::EC(InconsistentResolution::IgnoreInconsistent)
    };

    let h = butterfly::make_ecs(&folder.get_busfile(), mode);
    Ok(MyCUHistogram(h))  // conversion into a hashmap/dict
}


/// Per barcode, creates a CUHistogram: In that cell, how often do you see a UMI with n_reads
#[pyfunction]
pub (crate) fn make_ecs_across_cb(
    busfile: &str, matrix_file: &str, transcript_file: &str,
    // collapse_ec:bool,
) -> PyResult<HashMap<u64, MyCUHistogram>> {

    let folder = BusFolder::from_files(busfile, matrix_file, transcript_file);
    let mut h_read: HashMap<u64, MyCUHistogram> = HashMap::new();
    
    for (_cb, records) in folder.get_iterator().groupby_cb() {
        let mut h_read_cb = CUHistogram::new();
        for r in records {
            h_read_cb.add_counts(r.COUNT as usize, 1);
        }
        // let cu = CUHistogram::new(h_read_cb);
        // a bit annoying: we have to wrap it so we can return to python
        h_read.insert(_cb, MyCUHistogram(h_read_cb));
    }
    Ok(h_read)  // conversion into a hashmap/dict
}

/// Counts frequencies of barcodes per cell.
/// THIS IS NOT "PER BARCODE"
/// Returns two dictionaries (each key is a barcode):
/// - the frequency of a barcode in terms of UMIs
/// - the frequency of a barcode in terms of reads
#[pyfunction]
pub (crate) fn make_ecs_cb(
    busfile: &str, matrix_file: &str, transcript_file: &str,
    // collapse_ec:bool,
) -> PyResult<(MyCUHistogram, MyCUHistogram)> {

    let folder = BusFolder::from_files(busfile, matrix_file, transcript_file);
    
    // let mode = if collapse_ec {
    //     // println!("Making mapper");
    //     let mapper = folder.make_mapper(t2g_file);
    //     MappingMode::Gene(mapper, InconsistentResolution::IgnoreInconsistent)
    // } else {
    //     MappingMode::EC(InconsistentResolution::IgnoreInconsistent)
    // };
    // println!("Done Making mapper");

    // TODO work with the CUHistogram directly
    let mut cu_umi = CUHistogram::new();
    let mut cu_reads = CUHistogram::new();

    for (_cb, records) in folder.get_iterator().groupby_cb() {
        let nreads: u32 = records.iter().map(|r| r.COUNT).sum();
        cu_reads.add_counts(nreads as usize, 1);

        let numi = records.len();
        cu_umi.add_counts(numi, 1);
    }
    Ok(
        (MyCUHistogram(cu_umi), MyCUHistogram(cu_reads))  // conversion into a hashmap/dict
    ) 
}

/// Counts frequencies of barcodes.
/// Returns two dictionaries (each key is a barcode):
/// - the frequency of a barcode in terms of UMIs
/// - the frequency of a barcode in terms of reads
#[pyfunction]
pub (crate) fn make_ecs_ec(busfile: &str, matrix_file: &str, transcript_file: &str, t2g_file: &str, collapse_ec:bool, mapping_mode: &str)
 -> PyResult<(HashMap<String, MyCUHistogram>, MyCUHistogram, MyCUHistogram)> {
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
 -> (HashMap<String,MyCUHistogram>, MyCUHistogram, MyCUHistogram) {

    println!("Making mapper");
    let ecmapper = folder.make_mapper_transcript();
    println!("Done Making mapper");

    let mut h_read: HashMap<(TranscriptId, usize), usize> = HashMap::new();
    let mut h_read_multimapped = CUHistogram::new();
    let mut h_read_inconsistent=  CUHistogram::new();

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
                h_read_multimapped.add_counts(nreads as usize, 1);
            },
            MappingResultTranscript::Inconsistent => {
                let nreads: u32 = records.iter().map(|r| r.COUNT).sum();
                h_read_inconsistent.add_counts(nreads as usize, 1);
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
            h.add_counts(freq, count);

        } else {
            let mut inner = HashMap::new();
            inner.insert(freq, count);
            let h = CUHistogram::from(inner);
            final_read.insert(tid, h);

        }
    }

    // convert tid to actual trasncript name
    let mut final_with_string: HashMap<String, MyCUHistogram> = HashMap::new();
    for (tid, cu) in final_read {
        let ix = tid.0 as usize;
        final_with_string.insert(all_transcripts[ix].0.clone(), MyCUHistogram(cu));
    }

    (final_with_string, MyCUHistogram(h_read_inconsistent), MyCUHistogram(h_read_multimapped))
}

fn make_ecs_ec_genemapping(folder: BusFolder, t2g_file: &str, collapse_ec:bool
) -> (HashMap<String, MyCUHistogram>, MyCUHistogram, MyCUHistogram) {

    println!("Making mapper");
    let ecmapper = folder.make_mapper(t2g_file);
    println!("Done Making mapper");
    // let mapping_mode =  if collapse_ec{
    //      MappingMode::Gene(ecmapper, InconsistentResolution::IgnoreInconsistent)
    // } else {
    //     MappingMode::EC(InconsistentResolution::IgnoreInconsistent)
    // };

    // let mapping_mode = MappingMode::Gene(ecmapper, InconsistentResolution::IgnoreInconsistent);

    // let mut h_umi: HashMap<(GeneId, usize), usize> = HashMap::new();
    let mut h_read: HashMap<(GeneId, usize), usize> = HashMap::new();
    let mut h_read_multimapped = CUHistogram::new();
    let mut h_read_inconsistent=  CUHistogram::new();

    let all_genes = ecmapper.get_gene_list();

    // iterate over all CB/UMI, and if the records map to a single gene
    // add an entry of #reads per UMI into `h_read`
    println!("Iterating cb/umi");
    let spinner = get_spinner();
    let mut counter = 0;
    for (_cbumi, records) in folder.get_iterator().groupby_cbumi() {
    
        match find_consistent(&records, &ecmapper) {
            bustools::consistent_genes::MappingResult::SingleGene(geneid) => {
                let nreads: u32 = records.iter().map(|r| r.COUNT).sum();
                let v = h_read.entry((geneid, nreads as usize)).or_insert(0);
                *v+=1;
            },
            bustools::consistent_genes::MappingResult::Multimapped(_glist) => {
                let nreads: u32 = records.iter().map(|r| r.COUNT).sum();
                h_read_multimapped.add_counts(nreads as usize, 1);

            },
            bustools::consistent_genes::MappingResult::Inconsistent => {
                let nreads: u32 = records.iter().map(|r| r.COUNT).sum();
                h_read_inconsistent.add_counts(nreads as usize, 1);
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
            h.add_counts(freq, count);

        } else {
            let mut inner = HashMap::new();
            inner.insert(freq, count);
            let h = CUHistogram::from(inner);
            final_read.insert(geneid, h);

        }
        // let h= final_read.entry(geneid).or_insert()
    }
    // convert gene id to actual gene name
    let mut final_with_string: HashMap<String, MyCUHistogram> = HashMap::new();
    for (geneid, cu) in final_read {
        let ix = geneid.0 as usize;
        final_with_string.insert(all_genes[ix].0.clone(), MyCUHistogram(cu));
    }

    (final_with_string, MyCUHistogram(h_read_inconsistent), MyCUHistogram(h_read_multimapped))
}


#[pyfunction]
pub (crate) fn estimate_tgt<'py>(py: Python<'py>, 
    busfile: &str, matrix_file: &str, transcript_file: &str, t2g_file: &str,
    // collapse_ec:bool,
) -> PyResult<(pyo3_polars::PyDataFrame, Bound<'py, PyDict>, Bound<'py, PyDict>)> {

    let folder = BusFolder::from_files(busfile, matrix_file, transcript_file);
    let ecmapper = folder.make_mapper(t2g_file);
    let (df, multi_counter, inconsistent_counter) = _estimate_tgt(busfile, MappingMode::Gene(ecmapper, InconsistentResolution::IgnoreInconsistent));
    
    // multi_counter.to_object();
    let multi_counter_python = PyDict::new_bound(py);
    for (k, v) in multi_counter {
        let genes = k.into_iter().map(|g| g.0).collect_vec();
        let s = PyTuple::new_bound(py, genes.iter());
        multi_counter_python.set_item(s, v).unwrap();
    }

    let inconsistent_counter_python = PyDict::new_bound(py);
    for (k, v) in inconsistent_counter {
        let ecs = k.into_iter().map(|g| g.0).collect_vec();
        let s = PyTuple::new_bound(py, ecs.iter());
        inconsistent_counter_python.set_item(s, v).unwrap();
    }

    Ok(
        (pyo3_polars::PyDataFrame(df), multi_counter_python, inconsistent_counter_python)
    )
}

/// Transcript per transcript
/// 
/// Per CB/UMI, how often do we run into inconsitennt gene mappings.
/// This could be due to chimeras: The CB/UMI swapped to another mRNA and hence it's
pub fn _estimate_tgt(busfile: &str, mapping_mode: MappingMode) -> (DataFrame, HashMap<BTreeSet<Genename>, usize>, HashMap<BTreeSet<EC>, usize>) {

    let mut cb_counter_pure: HashMap<u64, usize> = HashMap::new();
    let mut cb_counter_inconsistent: HashMap<u64, usize> = HashMap::new();
    let mut cb_counter_multi: HashMap<u64, usize> = HashMap::new();
    let mut cbs_set = HashSet::new();

    // how often do we see a particular set of genes multimapped (i.e. we couldnt decided which one it is)
    let mut multimap_counter: HashMap<BTreeSet<Genename>, usize> = HashMap::new();

    // cant really report the inconsistent genes (would require solving a LinearProgram), just return the ECs that are inconsistent for new
    let mut inconsistent_counter: HashMap<BTreeSet<EC>, usize> = HashMap::new();


    let reader = BusReader::new(busfile);

    let mut pure = 0;
    let mut multimapped = 0;
    let mut inconsistent = 0;
    let mut total = 0;

    for ((_cb, _umi), recordlist) in reader.groupby_cbumi() {
        total += 1;
        cbs_set.insert(_cb);
        match &mapping_mode {
            // check if we can uniquely match those read to the same gene
            // if not its either multimapped or inconsistent (could be a CB/UMI collision)            
            MappingMode::Gene(ecmapper, _resolution_mode) => {
                match find_consistent(&recordlist, ecmapper) {
                    MappingResult::SingleGene(_) => {
                        // increment our histogram
                        pure += 1;
                        let v = cb_counter_pure.entry(_cb).or_insert(0);
                        *v+=1;
                    }
                    MappingResult::Multimapped(_genes) => {
                        multimapped += 1;            
                        let v = cb_counter_multi.entry(_cb).or_insert(0);
                        *v+=1;

                        if _genes.len() <= 1 {
                            println!("{recordlist:?}");
                            println!("{_genes:?}" );
                            // panic!("BLA")
                        }
                       
                        let genenames: BTreeSet<_> = _genes.into_iter().map(|gene_id|  ecmapper.resolve_gene_id(gene_id)).collect();
                        let v = multimap_counter.entry(genenames).or_insert(0);
                        *v += 1;
                    }
                    // inconsistent, i.e mapping to two distinct genes
                    // the reasonable thin
                    MappingResult::Inconsistent => {
                        inconsistent += 1;
                        let v = cb_counter_inconsistent.entry(_cb).or_insert(0);
                        *v+=1;    


                        let ecs: BTreeSet<EC> = recordlist.iter().map(|r| EC(r.EC)).collect();
                        let v = inconsistent_counter.entry(ecs).or_insert(0);
                        *v+=1;

                    },
                }
            },
            MappingMode::EC(_) => todo!(),
            MappingMode::Transcript(_, _) => todo!(),
        }
    }

    println!(
        "Total CB-UMI {total}, pure {pure} ({}%), Multimapped {multimapped} ({}%), Discarded/Inconsistent {inconsistent} ({}%)",
        100.0*(pure as f32) / (total as f32),
        100.0*(multimapped as f32) / (total as f32),
        100.0*(inconsistent as f32) / (total as f32)
    );


    let cbs = cbs_set.into_iter().collect_vec();
    let series_pure = cbs.iter().map(|c| *cb_counter_pure.get(c).unwrap_or(&0) as u64).collect_vec(); 
    let series_incons = cbs.iter().map(|c| *cb_counter_inconsistent.get(c).unwrap_or(&0) as u64).collect_vec(); 
    let series_multi = cbs.iter().map(|c| *cb_counter_multi.get(c).unwrap_or(&0) as u64).collect_vec(); 

    let df: DataFrame = df!(
        "CB" => &cbs,
        "pure" => series_pure,
        "inconsistent" => &series_incons,
        "multi" => &series_multi,
    ).unwrap();
    (df, multimap_counter, inconsistent_counter)
}


#[cfg(test)]
mod testing {
    use super::*;

    #[test]
    /// ensure that the CUhistogram remains unchanged
    fn test_insta_make_ecs(){
        let bus_file = "/home/michi/bus_testing/bus_output_short/output.corrected.sort.bus";
        let matrix_file = "/home/michi/bus_testing/bus_output_short/matrix.ec";
        let transcript_file = "/home/michi/bus_testing/bus_output_short/transcripts.txt";
        let t2g_file = "/home/michi/bus_testing/transcripts_to_genes.txt";

        let cu = make_ecs(bus_file, matrix_file, transcript_file, t2g_file, false).unwrap();
        insta::with_settings!({sort_maps => true}, {
            insta::assert_yaml_snapshot!(cu.0.get_histogram())
        });
    }

    #[test]
    fn test_insta_make_ecs_across_cb(){
        let bus_file = "/home/michi/bus_testing/bus_output_short/output.corrected.sort.bus";
        let matrix_file = "/home/michi/bus_testing/bus_output_short/matrix.ec";
        let transcript_file = "/home/michi/bus_testing/bus_output_short/transcripts.txt";

        let cu = make_ecs_across_cb(bus_file, matrix_file, transcript_file).unwrap();

        // need to convert CUHist->hashmap
        let converted: HashMap<u64, HashMap<usize, usize>> =cu.into_iter().map(|(cb, cu)| (cb, cu.0.get_histogram()) ).collect();
        insta::with_settings!({sort_maps => true}, {
            insta::assert_yaml_snapshot!(converted)
        });
    }    

    #[test]
    fn test_insta_make_ecs_cb(){
        let bus_file = "/home/michi/bus_testing/bus_output_short/output.corrected.sort.bus";
        let matrix_file = "/home/michi/bus_testing/bus_output_short/matrix.ec";
        let transcript_file = "/home/michi/bus_testing/bus_output_short/transcripts.txt";

        let (cu_umi, cu_read) = make_ecs_cb(bus_file, matrix_file, transcript_file).unwrap();

        println!("{:?}", cu_umi.0);
        // need to convert CUHist->hashmap
        insta::with_settings!({sort_maps => true}, {
            insta::assert_yaml_snapshot!((cu_umi.0.get_histogram(), cu_read.0.get_histogram()))
        });
    }    

    #[test]
    fn test_insta_make_ecs_ec(){
        let bus_file = "/home/michi/bus_testing/bus_output_short/output.corrected.sort.bus";
        let matrix_file = "/home/michi/bus_testing/bus_output_short/matrix.ec";
        let transcript_file = "/home/michi/bus_testing/bus_output_short/transcripts.txt";
        let t2g_file = "/home/michi/bus_testing/transcripts_to_genes.txt";

        // gene mapping
        let (cu_per_ec, cu_incon, cu_multi) = make_ecs_ec(bus_file, matrix_file, transcript_file, t2g_file, false, "gene").unwrap();
        let converted: HashMap<String, HashMap<usize, usize>> =cu_per_ec.into_iter().map(|(cb, cu)| (cb, cu.0.get_histogram()) ).collect();
        insta::with_settings!({sort_maps => true}, {
            insta::assert_yaml_snapshot!(converted);
            insta::assert_yaml_snapshot!(cu_incon.0.get_histogram());
            insta::assert_yaml_snapshot!(cu_multi.0.get_histogram());
        });

        // transcript mapping
        let (cu_per_ec, cu_incon, cu_multi) = make_ecs_ec(bus_file, matrix_file, transcript_file, t2g_file, false, "transcript").unwrap();
        let converted: HashMap<String, HashMap<usize, usize>> =cu_per_ec.into_iter().map(|(cb, cu)| (cb, cu.0.get_histogram()) ).collect();
        insta::with_settings!({sort_maps => true}, {
            insta::assert_yaml_snapshot!(converted);
            insta::assert_yaml_snapshot!(cu_incon.0.get_histogram());
            insta::assert_yaml_snapshot!(cu_multi.0.get_histogram());
        });
    }   
}