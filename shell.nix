{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell {
  buildInputs = with pkgs; [
    bashInteractive
    rEnv
    rstudioEnv
    gsl
    git
    rPackages.devtools
    rPackages.rmarkdown
    rPackages.ggdist
    rPackages.hypergeo
    pandoc
    # ... any other dependencies
  ];

  # You can set environment variables as needed
  # environment variables that should be set when the shell is entered
  shellHook = ''
    echo "entering rEnv rstudioEnv and gsl evnironment"
    mkdir -p "$(pwd)/_libs"
    export R_LIBS_USER="$(pwd)/_libs"
    #  export SOME_VAR=some_value
  '';
}
