test.desc.cats <- function() {
    cats <- get.desc.categories()
    print(cats)
    checkEquals(5, length(cats))
}

test.desc.names <- function() {
    cats <- get.desc.categories()
    for (acat in cats) {
        dnames <- get.desc.names(acat)
        checkTrue(length(dnames) > 0)
    }
}

test.desc.calc <- function() {
    dnames <- get.desc.names("topological")
    mols <- parse.smiles("c1ccccc1CCC")
    dvals <- eval.desc(mols, dnames[1])
    checkTrue(dvals[1,1] == 1)
}

test.decs.parameter.passing <- function() {
    # Get a descriptor with parameters: 'LongestAliphaticChain'
    dname <- get.desc.names()[20]

    # This molecules races an StackOverflow exception, if we do not set the
    # parameter: 'checkRingSystems' = TRUE, for 'LongestAliphaticChain'
    # Default for it is FALSE
    mol <- parse.smiles("c1ccccc1C")[[1]]
    dvals <- eval.desc(mol, dname, desc.params=list(checkRingSystem=TRUE))
    checkEqualsNumeric(dvals[1,1], 0)

    mol <- parse.smiles("c1ccccc1C")[[1]]
    dvals <- eval.desc(mol, dname, desc.params=list(checkRingSystem=FALSE))
    checkTrue(is.na(dvals[1,1]))

    mol <- parse.smiles("c1ccccc1C")[[1]]
    dvals <- eval.desc(mol, dname)
    checkTrue(is.na(dvals[1,1]))

    # Test empty parameter list
    mol <- parse.smiles("C=C(CCC1CC1C(C)C(C)C)C(C)CC2CCCC2")[[1]]
    dvals <- eval.desc(mol, get.desc.names()[1], desc.params=list())
    checkEqualsNumeric(dvals[1,1], 0.9)

    mol <- parse.smiles("c1ccccc1C")[[1]]
    dvals <- eval.desc(mol, get.desc.names()[20], desc.params=list())
    checkTrue(is.na(dvals[1,1]))

    # Test wrong parameter list
    mol <- parse.smiles("C=C(CCC1CC1C(C)C(C)C)C(C)CC2CCCC2")[[1]]
    checkException(eval.desc(mol, get.desc.names()[41],
                             desc.params=list(nhigh=as.integer(1))))  # wrong size
    checkException(eval.desc(mol, get.desc.names()[41],
                             desc.params=list(nhigh=as.integer(1),
                                              nlow=1, # <-- wrong type
                                              checkAromaticity=TRUE)))

    # Test that order is internally corrected
    dvals1 <- eval.desc(mol, get.desc.names()[41],
                        desc.params=list(nlow=as.integer(2),
                                         nhigh=as.integer(1),
                                         checkAromaticity=FALSE))
    dvals2 <- eval.desc(mol, get.desc.names()[41],
                        desc.params=list(checkAromaticity=FALSE,
                                         nlow=as.integer(2),
                                         nhigh=as.integer(1)))
    dvals3 <- eval.desc(mol, get.desc.names()[41],
                        desc.params=list(checkAromaticity=FALSE,
                                         nhigh=as.integer(1),
                                         nlow=as.integer(2)))

    checkEquals(dvals1, dvals2)
    checkEquals(dvals1, dvals3)
}

test.desc.lac <- function() {
    dname <- get.desc.names()[20]  # LongestAliphaticChain

    mol <- parse.smiles("CCCCc1ccccc1")[[1]]
    dvals <- eval.desc(mol, dname, desc.params=list(checkRingSystem=TRUE))
    checkEqualsNumeric(dvals[1,1], 4)

    mol <- parse.smiles("C=C(CCC1CC1C(C)C(C)C)C(C)CC2CCCC2")[[1]]
    dvals <- eval.desc(mol, dname, desc.params=list(checkRingSystem=TRUE))
    checkEqualsNumeric(dvals[1,1], 5)

    mol <- parse.smiles("CCO")[[1]]
    dvals <- eval.desc(mol, dname, desc.params=list(checkRingSystem=TRUE))
    checkEqualsNumeric(dvals[1,1], 2)

    # Note here: The 'LongestAliphaticChain' descriptor modifies the atom
    # container if 'checkRingSystem'= TRUE. So better one clones the container
    # before calculating the descriptors with modified parameters.
}