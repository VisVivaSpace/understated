//! Example: Query ephemeris data from SPK files
//!
//! Demonstrates using the Ephemeris API to query planetary positions and
//! velocities at arbitrary epochs with automatic body chaining.
//!
//! Run with: cargo run --example query_ephemeris --features test-data -- test_data/de440s.bsp

use std::env;
use understated::{Ephemeris, EpochTDB, NaifId, Rectangular};
use understated::time::format_iso8601;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        println!("=== SPK Ephemeris Query Example ===\n");
        println!("Usage: {} <spk_file.bsp> [epoch]", args[0]);
        println!();
        println!("Arguments:");
        println!("  spk_file.bsp  - SPK ephemeris file");
        println!("  epoch         - Optional: time string (default: J2000.0)");
        println!();
        println!("Examples:");
        println!("  {} de440s.bsp", args[0]);
        println!("  {} de440s.bsp \"2020-06-15T12:00:00\"", args[0]);
        println!();
        println!("API Overview:\n");
        println!("  let eph = Ephemeris::load(path)?;");
        println!("  let state = eph.state_of(target, epoch, center)?;");
        println!();
        println!("  // State contains:");
        println!("  //   target:   NaifId      // Body this state describes");
        println!("  //   center:   NaifId      // Origin of coordinates");
        println!("  //   frame:    i32         // Reference frame (e.g., 1 = J2000)");
        println!("  //   position: [f64; 3]    // km");
        println!("  //   velocity: [f64; 3]    // km/s");

        return Ok(());
    }

    let spk_path = &args[1];
    let epoch = if args.len() > 2 {
        EpochTDB::parse(&args[2])?
    } else {
        EpochTDB(0.0) // J2000.0
    };

    println!("=== SPK Ephemeris Query Example ===\n");
    println!("SPK file: {}", spk_path);
    println!(
        "Epoch: {} (TDB = {:.3} s)\n",
        format_iso8601(epoch.0),
        epoch.0
    );

    // Load the kernel
    let eph = Ephemeris::load(spk_path)?;

    // List available bodies
    let bodies = eph.spk_bodies();
    println!("Available bodies: {:?}\n", bodies);

    // Query known bodies relative to SSB
    let queries: &[(NaifId, &str)] = &[
        (NaifId::EARTH, "Earth"),
        (NaifId::MOON, "Moon"),
        (NaifId::MARS, "Mars"),
        (NaifId::SUN, "Sun"),
    ];

    println!("States relative to Solar System Barycenter:\n");

    for (target, name) in queries {
        match eph.state_of(*target, epoch, NaifId::SSB) {
            Ok(state) => {
                let pos = Rectangular(state.position);
                let lat = pos.to_latitudinal();

                println!(
                    "  {} (NAIF {}) -> SSB (frame {}):",
                    name, state.target, state.frame
                );
                println!(
                    "    Position: [{:>15.3}, {:>15.3}, {:>15.3}] km",
                    state.position[0], state.position[1], state.position[2]
                );
                println!(
                    "    Velocity: [{:>15.6}, {:>15.6}, {:>15.6}] km/s",
                    state.velocity[0], state.velocity[1], state.velocity[2]
                );
                println!(
                    "    Distance: {:.3} km ({:.6} AU)",
                    state.distance(),
                    state.distance() / 149597870.7
                );
                println!("    Speed: {:.6} km/s", state.speed());
                println!(
                    "    Lat/Lon: {:.3} deg / {:.3} deg",
                    lat.latitude_deg(),
                    lat.longitude_deg()
                );
                println!();
            }
            Err(e) => {
                println!("  {}: not available ({:?})\n", name, e);
            }
        }
    }

    // Demonstrate body chaining: Moon relative to Earth
    println!("Body chaining demonstration:");
    println!("  (Moon relative to Earth, requires chaining through EMB)\n");

    match eph.state_of(NaifId::MOON, epoch, NaifId::EARTH) {
        Ok(state) => {
            println!(
                "  Moon-Earth distance: {:.3} km",
                state.distance()
            );
            println!("  Moon-Earth speed: {:.6} km/s", state.speed());
        }
        Err(e) => {
            println!("  Could not query Moon relative to Earth: {:?}", e);
        }
    }

    Ok(())
}
