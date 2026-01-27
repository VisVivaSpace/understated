use std::fs::File;
use muad_dib::{DAFFile, DAFSegment};

fn main() {
    let file = File::open("test_data/test.bc").expect("Failed to open test.bc");
    let daf = DAFFile::from_file(file).expect("Failed to parse CK DAF");

    println!("=== test.bc CK Segments ===\n");
    
    let mut count = 0;
    for (i, seg_result) in daf.enumerate() {
        if let Ok(segment) = seg_result {
            if let DAFSegment::CK(ck) = segment {
                count += 1;
                println!("Segment {}:", i);
                println!("  Instrument Code: {}", ck.instrument_code);
                println!("  CK Type: {}", ck.ck_type);
                println!("  Frame Code: {}", ck.frame_code);
                println!("  Initial SCLK: {}", ck.initial_sclk);
                println!("  Final SCLK: {}", ck.final_sclk);
                println!("  Has Angular Rates: {}", ck.rates);
                println!("  Data Size: {}", ck.data.len());
                println!();
            }
        }
    }
    
    println!("Total CK segments: {}", count);
}
