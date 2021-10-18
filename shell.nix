{}:
let
    pkgs = import <nixpkgs> {};
    mgc = (import ./default.nix) {};
in pkgs.mkShell {
    name = "manif-geom-cpp-env";
    src = null;
    shellHook = ''
        export PATH=$PATH:${mgc}/bin
    '';
}