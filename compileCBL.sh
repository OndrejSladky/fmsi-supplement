cd CBL
git pull

echo "======== compiling CBL for k=15 =============="
RUSTFLAGS="-C target-cpu=native" K=15 PREFIX_BITS=28 cargo +nightly build --release --examples --target-dir target.k_15

echo "======== compiling CBL for k=23 =============="
RUSTFLAGS="-C target-cpu=native" K=23 PREFIX_BITS=28 cargo +nightly build --release --examples --target-dir target.k_23

echo "======== compiling CBL for k=31 =============="
RUSTFLAGS="-C target-cpu=native" K=31 PREFIX_BITS=28 cargo +nightly build --release --examples --target-dir target.k_31