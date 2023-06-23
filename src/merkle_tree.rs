use rs_merkle::{MerkleTree, MerkleProof, algorithms::Sha256, Hasher};
struct MerkleTreeProver {
    leaf_values: Vec<&'static str>,
    leaves: Vec<[u8; 32]>,
    merkle_tree: MerkleTree<Sha256>
}

struct MerkleTreeVerifier { 
    merkle_root: [u8; 32],
    leave_number: usize
}

impl MerkleTreeProver {
    fn new(leaf_values: Vec<&'static str>) -> Self {
        let leaves: Vec<[u8; 32]> = leaf_values.iter()
            .map(|x| Sha256::hash(x.as_bytes())).collect();
        let merkle_tree = MerkleTree::<Sha256>::from_leaves(&leaves);
        Self {
            leaf_values,
            leaves,
            merkle_tree
        }
    }

    fn commit(&self) -> [u8; 32] {
        self.merkle_tree.root().unwrap()
    }

    fn open(&self, leaf_indices: &Vec<usize>) -> Vec<u8> {
        self.merkle_tree.proof(leaf_indices).to_bytes()
    }
}

impl MerkleTreeVerifier {
    fn new(leave_number: usize, merkle_root: &[u8; 32]) -> Self {
        Self { 
            leave_number, 
            merkle_root: *merkle_root 
        }
    }

    fn verify(&self, proof_bytes: Vec<u8>, indices: Vec<usize>, leaves: &Vec<&str>) {
        let proof = MerkleProof::<Sha256>::try_from(proof_bytes).unwrap();
        let leaves_to_prove: Vec<[u8; 32]> = leaves.iter()
            .map(|x| Sha256::hash(x.as_bytes())).collect();
        assert!(proof.verify(self.merkle_root, &indices, &leaves_to_prove, self.leave_number));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn commit_and_open() {
        let leaf_values = vec!["a", "b", "c", "d", "e", "f"];
        let leave_number = leaf_values.len();
        let prover = MerkleTreeProver::new(leaf_values);
        let root = prover.commit();
        let verifier = MerkleTreeVerifier::new(leave_number, &root);
        let leaf_indices = vec![2, 3];
        let proof_bytes = prover.open(&leaf_indices);
        verifier.verify(proof_bytes, leaf_indices, &vec!["c", "d"]);
    }
}
