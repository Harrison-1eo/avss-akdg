use rs_merkle::{MerkleTree, MerkleProof, Hasher};


#[derive(Clone)]
pub struct Blake3Algorithm {}

impl Hasher for Blake3Algorithm {
    type Hash = [u8; 32];

    fn hash(data: &[u8]) -> [u8; 32] {
        blake3::hash(data).into()
    }
}

struct MerkleTreeProver {
    merkle_tree: MerkleTree<Blake3Algorithm>
}

struct MerkleTreeVerifier { 
    merkle_root: [u8; 32],
    leave_number: usize
}

impl MerkleTreeProver {
    fn new(leaf_values: Vec<Vec<u8>>) -> Self {
        let leaves: Vec<[u8; 32]> = leaf_values.iter()
            .map(|x| Blake3Algorithm::hash(x)).collect();
        let merkle_tree = MerkleTree::<Blake3Algorithm>::from_leaves(&leaves);
        Self {
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
            merkle_root: merkle_root.clone()
        }
    }

    fn verify(&self, proof_bytes: Vec<u8>, indices: Vec<usize>, leaves: &Vec<Vec<u8>>) -> bool {
        let proof = MerkleProof::<Blake3Algorithm>::try_from(proof_bytes).unwrap();
        let leaves_to_prove: Vec<[u8; 32]> = leaves.iter()
            .map(|x| Blake3Algorithm::hash(x)).collect();
        proof.verify(self.merkle_root, &indices, &leaves_to_prove, self.leave_number)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::field::{mersenne61_ext::Mersenne61Ext, Field, as_bytes_vec};

    #[test]
    fn commit_and_open() {
        let leaf_values = vec![
            as_bytes_vec(&[Mersenne61Ext::from_int(1), Mersenne61Ext::from_int(2)]), 
            as_bytes_vec(&[Mersenne61Ext::from_int(3), Mersenne61Ext::from_int(4)]), 
            as_bytes_vec(&[Mersenne61Ext::from_int(5), Mersenne61Ext::from_int(6)]), 
            as_bytes_vec(&[Mersenne61Ext::from_int(7), Mersenne61Ext::from_int(8)]), 
            as_bytes_vec(&[Mersenne61Ext::from_int(9), Mersenne61Ext::from_int(10)]), 
            as_bytes_vec(&[Mersenne61Ext::from_int(11), Mersenne61Ext::from_int(12)]), 
            as_bytes_vec(&[Mersenne61Ext::from_int(13), Mersenne61Ext::from_int(14)]), 
        ];
        let leave_number = leaf_values.len();
        let prover = MerkleTreeProver::new(leaf_values);
        let root = prover.commit();
        let verifier = MerkleTreeVerifier::new(leave_number, &root);
        let leaf_indices = vec![2, 3];
        let proof_bytes = prover.open(&leaf_indices);
        let open_values = vec![
            as_bytes_vec(&[Mersenne61Ext::from_int(5), Mersenne61Ext::from_int(6)]),
            as_bytes_vec(&[Mersenne61Ext::from_int(7), Mersenne61Ext::from_int(8)])
        ];
        verifier.verify(proof_bytes, leaf_indices, &open_values);
    }

    #[test]
    fn blake3() {
        let hash_res = Blake3Algorithm::hash("data".as_bytes());
        let hex_string = hex::encode(hash_res);
        assert_eq!("28a249c2e4d3a92bc0a16ed8f1b5cf83ca20415ee12e502b096624902bbc97bd", hex_string);
    }
}