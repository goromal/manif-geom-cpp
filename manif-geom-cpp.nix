{ stdenv
, cleanSource
, cmake
, clang
, git
, eigen
, boost
}:
stdenv.mkDerivation {
    name = "manif-geom-cpp";
    version = "1.0.0";
    src = cleanSource ./.;
    nativeBuildInputs = [
        cmake
        clang
        git
    ];
    buildInputs = [
        eigen
        boost
    ];
}