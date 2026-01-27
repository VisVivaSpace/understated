fn main() {
    if cfg!(feature = "cspice") {
        println!("cargo:rustc-link-search=native=/Users/nstrange/cspice/lib");
    }
}
