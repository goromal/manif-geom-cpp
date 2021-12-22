{ pkgs ? import <nixpkgs> {} }:
import ./manif-geom-cpp.nix {
    stdenv = pkgs.clangStdenv;
    cleanSource = pkgs.lib.cleanSource;
    cmake = pkgs.cmake;
    clang = pkgs.clang;
    git = pkgs.git;
    eigen = pkgs.eigen;
    boost = pkgs.boost;
}

