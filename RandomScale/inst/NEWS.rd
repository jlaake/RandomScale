\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}

\name{NEWS}
\title{RandomScale News}
\encoding{UTF-8}
\section{Changes in version 0.1.6 (0.1.5 skipped)(2013-9-18)}{
  \itemize{
    \item Added compute.Nhat and compute.Nhat.se
  }
}

\section{Changes in version 0.1.4 (2013-3-5)}{
  \itemize{
    \item Modified tpl files to use M_PI instead of PI
    \item Fixed code in fitadmb when likeihood="fixed" was incorrectly testing for sigma
  }
}

\section{Changes in version 0.1.3 (2012-12-21)}{
  \itemize{
    \item Modified tpl files to use mfexp and estimation phases; estimation phases are not used for g-likelihood where it caused problems
    \item Added fixed effect model to fitadmb; likelihood="fixed"
    \item Added nsteps argument for fixed effect fit; number of steps for adromb integration
    \item Distances are scaled (x/max(x)) for estimation but parameter and likelihood are adjusted to original scale
    \item In fitdata, added starting value computation for beta, bounds on parameters, added method argument with default L-BFGS-B, added debug argument
    \item In fitadmb, added starting value computations with creation of PIN file, added debug argument; added keep argument to keep current exe and tpl; set default verbose=TRUE
  }
}


\section{Changes in version 0.1.2 (2012-12-11)}{
  \itemize{
    \item in fitadmb, call clean_admb before run_admb to make sure no output files are left over that could mess up read_admb
  }
}

\section{Changes in version 0.1.1 (2012-12-05)}{
  \itemize{
    \item Patches to documentation
    \item Moved example DAT files from Data to Inst/Data
  }
}


\section{Changes in version 0.1.0 (2012-12-02)}{
  \itemize{
    \item Addition of mixed effects modeling through ADMB
    \item Added examples to help
  }
}

\section{Initial release version 0.0.1}{
  \itemize{
    \item Random effects half-normal scale from R and ADMB
  }
}

