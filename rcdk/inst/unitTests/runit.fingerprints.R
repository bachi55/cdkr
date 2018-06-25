test.fp <- function() {
    mol <- parse.smiles("CCCCC")[[1]]
    fp <- get.fingerprint(mol, type='maccs')
    checkTrue(length(fp@bits) > 0)
    fp <- get.fingerprint(mol, type='kr')
    checkTrue(length(fp@bits) > 0)
    fp <- get.fingerprint(mol, type='shortestpath')
    checkTrue(length(fp@bits) > 0)
}

test.fp.substructures <- function() {
    # User defined patterns
    smarts <- c("c1ccccc1", "[CX4H3][#6]", "[CX2]#[CX2]")
    mol <- parse.smiles("c1ccccc1CCC")[[1]]
    fp <- get.fingerprint(mol, type="substructure")
    
    checkEquals(length(fp), 3)
    checkEquals(length(fp@bits), 2)
    checkEquals(fp@bits[1], 1)
    checkEquals(fp@bits[2], 2)
    
    
    mol <- parse.smiles("C=C=C")[[1]]
    fp <- get.fingerprint(mol, type="substructure")
    
    checkEquals(length(fp), 3)
    checkEquals(length(fp@bits), 0)
}
