\name{eval.desc}
\alias{eval.desc}
\title{
  Evaluate a Molecular Descriptor
}
\description{
  The CDK implements a number of descriptors divided into three main
  groups - atomic, molecular and bond. This method evaluates the
  specified molecular  descriptor(s) for a molecule
}
\usage{
eval.desc(molecules, which.desc, verbose=FALSE)
}
\arguments{
  \item{molecules}{A single \code{IAtomContainer} object or a  \code{list} of references to CDK
  \code{IAtomContainer} objects}
  \item{which.desc}{The fully qualified class name of the descriptor
  to evaluate or a vector such names}
  \item{verbose}{If TRUE, progress will be written to the screen,
  otherwise the function performs silently}
  \item{desc.params}{Names list of parameters passed to the chosen descriptor. Currently this parameter can
  can only be used, if a single descriptor should be calculates, i.e. \code{length(which.desc)==1}.}
}
\examples{
smiles <- c('CCC', 'c1ccccc1', 'CC(=O)C')
mols <- sapply(smiles, parse.smiles)

dnames <- get.desc.names('constitutional')
descs <- eval.desc(mols, dnames, verbose=TRUE)

### Passing descriptor parameters to LongestAliphaticChainDescriptor
# For some molecules this descriptor would fail, if the atom container is not
# checked for rings:
mol <- parse.smiles("c1ccccc1C")[[1]]
dvals <- eval.desc(mol, get.desc.names()[20], desc.params=list(checkRingSystem=TRUE))

# PLEASE NOTE: The 'LongestAliphaticChain' descriptor modifies the atom
# container if 'checkRingSystem'= TRUE. So better one clones the container
# before calculating the descriptors with modified parameters. You can use, e.g.:
mol_copy <- .jcall(mol, "Lorg/openscience/cdk/interfaces/IAtomContainer;", "clone")
}
\value{
	A \code{data.frame} is returned. For a single molecule it will
	have one row, for multiple molecules it will have the
	number of rows equal to the number of molecules
}
\seealso{
\code{\link{get.desc.names}}
\code{\link{get.desc.categories}}
}
\keyword{programming}

\author{Rajarshi Guha (\email{rajarshi.guha@gmail.com})}
