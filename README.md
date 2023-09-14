# Benchmarking Froll and FRAVSS

## Overall

This repo implements a series of  FRI-based multilinear polynomial commitment scheme such as Gemini, Virgo and Froll, which are mainly in directory `gemini-fri/src`, `virgo/src` and `pcs/src` respectively.

This repo also implements FRI-based one to many univariate polynomial and binary polynomial commitment scheme, which are mainly in directory `vss/src` and `avss/src` respectively.

All protocols mentioned above use `util`, which implements merkle tree, finite field, polynomial and other utilities.

## Setup

Install `Rust`. See [Rust Installation](https://www.rust-lang.org/tools/install) for detail.

After installation, running `cargo --version` and `rustup --version` can see version of cargo and rustup.

Use Rust nightly toolchain by running `rustup default nightly`.

## Benchmarking All Schemes

Run `cargo bench`. It will benchmarking all schemes.

## Run Tests and Get Proof Size

Run `cargo test -- --nocapture`. It will also output proof sizes for protocols.



