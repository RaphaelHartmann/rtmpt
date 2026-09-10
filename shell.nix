let
  pkgs = import (fetchTarball "https://github.com/NixOS/nixpkgs/archive/8c3cede7ddc26bd659d2d383b5610efbd2c7a16e.tar.gz") {};
  lib = pkgs.lib;
in

with pkgs;

let
  rPkgs = with rPackages; [
    coda
    data_table
    loo
    Ryacas
    stringr
    truncnorm
    knitr
    rmarkdown
    devtools
    roxygen2
    testthat
    ggdist
    hypergeo
  ];
  sitePkgs = buildEnv {
    name = "rtmpt-r-site";
    paths = lib.closePropagation rPkgs;
    pathsToLink = [ "/library" ];
    ignoreCollisions = true;
  };

  # Build rr from master branch for Arrow Lake (Core Ultra 5 225U) support.
  # rr 5.9.0 does not recognize CPU ID 0xb0650; upstream master already has the fix.
  R = pkgs.R.overrideAttrs (old: {
    NIX_CFLAGS_COMPILE = (old.NIX_CFLAGS_COMPILE or "") + " -g -O0";
    NIX_LDFLAGS = (old.NIX_LDFLAGS or "") + " -g";
    dontStrip = true;
    doCheck = false;
  });

  rr = pkgs.rr.overrideAttrs (old: {
    name = "rr-git-master";
    src = builtins.fetchTarball {
      url = "https://github.com/rr-debugger/rr/archive/master.tar.gz";
    };
    patches = []; # master already has all upstream fixes
    # Dependencies from the nixpkgs rr derivation are reused as-is.
  });
in

mkShell {
  buildInputs = [
    bashInteractive
    rstudio
    R
    gsl
    gsl.dev
    git
    pandoc
    which
    gfortran
    cmake
    pkg-config
    gdb
    gdbgui
    radian
    valgrind
    rr
  ] ++ rPkgs;

  shellHook = ''
    export R_LIBS_USER="$(pwd)/_libs"
    mkdir -p "$R_LIBS_USER"
    export R_LIBS_SITE="${sitePkgs}/library"
    if [ -n "$LD_LIBRARY_PATH" ]; then
      export LD_LIBRARY_PATH="${gsl}/lib:$LD_LIBRARY_PATH"
    else
      export LD_LIBRARY_PATH="${gsl}/lib"
    fi
    export GSL_CFLAGS="$(gsl-config --cflags)"
    export GSL_LIBS="$(gsl-config --libs)"
    mkdir -p "$R_LIBS_USER/Makevars"
    cat > "$R_LIBS_USER/Makevars/Makevars" <<'RMAKEVARS'
GSL_CFLAGS = $(shell gsl-config --cflags)
GSL_LIBS = $(shell gsl-config --libs)
PKG_CFLAGS = $(GSL_CFLAGS)
PKG_CXXFLAGS = $(GSL_CFLAGS)
PKG_LIBS = $(GSL_LIBS)
CXX_STD = CXX11
CXXFLAGS += -ggdb1 -O0 -fno-omit-frame-pointer
CFLAGS   += -ggdb1 -O0 -fno-omit-frame-pointer
LDFLAGS  += -ggdb1
RMAKEVARS
    cat > "$R_LIBS_USER/Makevars/Makevars.asan" <<'RMAKEVARS_ASAN'
GSL_CFLAGS = $(shell gsl-config --cflags)
GSL_LIBS = $(shell gsl-config --libs)
PKG_CFLAGS = $(GSL_CFLAGS) -fsanitize=address -fno-omit-frame-pointer
PKG_CXXFLAGS = $(GSL_CFLAGS) -fsanitize=address -fno-omit-frame-pointer
PKG_LIBS = $(GSL_LIBS) -fsanitize=address
CXX_STD = CXX11
CXXFLAGS += -ggdb1 -O0 -fno-omit-frame-pointer
CFLAGS   += -ggdb1 -O0 -fno-omit-frame-pointer
LDFLAGS  += -ggdb1
RMAKEVARS_ASAN
    export R_MAKEVARS_USER="$R_LIBS_USER/Makevars/Makevars"
    echo
    printf '\033[1;36m╔═══════════════════════════════════════════════╗\033[0m\n'
    printf '\033[1;36m║    \033[1;33mrtmpt R development shell\033[1;36m                 ║\033[0m\n'
    printf '\033[1;36m╠═══════════════════════════════════════════════╣\033[0m\n'
    printf '\033[1;36m║  \033[0mDebug:  \033[1;33mgdb\033[0m, \033[1;33mvalgrind\033[0m, \033[1;33mrr\033[0m installed        \033[1;36m║\033[0m\n'
    printf '\033[1;36m║  \033[0mFlags:  \033[1;33m-O1 -ggdb -fno-omit-frame-pointer\033[0m   \033[1;36m║\033[0m\n'
    printf '\033[1;36m║  \033[0mASan:   \033[1;33mMakevars.asan\033[0m available           \033[1;36m║\033[0m\n'
    printf '\033[1;36m║  \033[0mUsage:  \033[1;33mexport R_MAKEVARS_USER="\$R_LIBS_USER/Makevars/Makevars.asan"\033[0m \033[1;36m║\033[0m\n'
    printf '\033[1;36m║  \033[0mThen    \033[1;33mR CMD INSTALL rtmpt\033[0m                 \033[1;36m║\033[0m\n'
    printf '\033[1;36m║  \033[0mRun     \033[1;33mR -d rr\033[0m or \033[1;33mR -d gdb\033[0m              \033[1;36m║\033[0m\n'
    printf '\033[1;36m║  \033[0mRecord: \033[1;33mrr record -- R -e "library(rtmpt); ..."\033[0m \033[1;36m║\033[0m\n'
    printf '\033[1;36m╚═══════════════════════════════════════════════╝\033[0m\n'
    echo
    PS1='\033[1;36m[rtmpt]\033[0m \033[1;33m\u@\h\033[0m:\033[1;34m\w\033[0m\$ '
  '';
}
