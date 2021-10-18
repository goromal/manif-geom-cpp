{}:
let
    pkgs = import <nixpkgs> {};
in pkgs.clangStdenv.mkDerivation {
    name = "manif-geom-cpp";
    version = "0.0.0";
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

