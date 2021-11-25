{ pkgs ? import <nixpkgs> {} }:
pkgs.clangStdenv.mkDerivation {
    name = "manif-geom-cpp";
    version = "1.0.0";
    src = pkgs.lib.cleanSource ./.;
    
    nativeBuildInputs = [
        pkgs.cmake
        pkgs.clang
        pkgs.git
    ];
    
    buildInputs = [
        pkgs.eigen
        pkgs.boost
    ];
}

