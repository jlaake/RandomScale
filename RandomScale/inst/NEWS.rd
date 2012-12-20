\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}

\name{NEWS}
\title{RandomScale News}
\encoding{UTF-8}

\section{Changes in version 0.1.3 (2012-12-19)}{
  \itemize{
    \item Modified tpl files to use mfexp and estimation phases; estimation phases are not used for g-likelihood where it caused problems
    \item Added fixed effect model to fitadmb; likelihood="fixed"
    \item Added nsteps argument for fixed effect fit; number of steps for adromb integration
    \item Distances are scaled (x/max(x)) for estimation but parameter and likelihood are adjusted to original scale
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

