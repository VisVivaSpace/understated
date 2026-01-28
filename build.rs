fn main() {
    if cfg!(feature = "cspice") {
        let cspice_dir = std::env::var("CSPICE_DIR")
            .expect("CSPICE_DIR environment variable must be set when building with the `cspice` feature");
        println!("cargo:rustc-link-search=native={}/lib", cspice_dir);
    }
}
