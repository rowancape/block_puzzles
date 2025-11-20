use ahash::AHashSet;

use crate::common::{Block, Coordinate};

pub struct Rotator {
    pub block: Block,
    pub axis_rot_state: [u8; 3],
    pub previous_rotations: AHashSet<Block>
}

impl Rotator {
    pub fn new(block: Block) -> Self {
        Self {
            block,
            axis_rot_state: [0, 0, 0],
            previous_rotations: AHashSet::new()
        }
    }
 
    pub fn rotate_by_i(&mut self, i: usize) {
        if !self.previous_rotations.contains(&self.block) {
            self.previous_rotations.insert(self.block.clone());
        }

        self.axis_rot_state[i] = (self.axis_rot_state[i] + 1) % 4;
        let max = self.block.dimensions.clone();

        match i {
            0 => { self.rotate_x(&max); }
            1 => { self.rotate_y(&max); }
            2 => { self.rotate_z(&max); }
            _ => { panic!("Invalid rotation index in Rotator function 'rotate_by_i()'!") }
        }
    }

    fn rotate_x(&mut self, max: &Coordinate) {
        for coord in &mut self.block.body {
            let y = coord.y;

            coord.y = coord.z;
            coord.z = (max.y - 1) - y;
        }

        let temp = self.block.dimensions.y;
        self.block.dimensions.y = self.block.dimensions.z;
        self.block.dimensions.z = temp; 
    }

    fn rotate_y(&mut self, max: &Coordinate) {
        for coord in &mut self.block.body {
            let x = coord.x;

            coord.x = coord.z;
            coord.z = max.x - x - 1;
        }

        let temp = self.block.dimensions.x;
        self.block.dimensions.x = self.block.dimensions.z;
        self.block.dimensions.z = temp; 
    }

    fn rotate_z(&mut self, max: &Coordinate) {
        for coord in &mut self.block.body {
            let x = coord.x;

            coord.x = coord.y;
            coord.y = max.x - x - 1;
        }

        let temp = self.block.dimensions.y;
        self.block.dimensions.y = self.block.dimensions.x;
        self.block.dimensions.x = temp; 
    }

    pub fn is_current_rot_unique(&self) -> bool {
        !self.previous_rotations.contains(&self.block)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rotation_dimensions() {
        let block = Block::new(vec![vec![
            vec![1, 1, 1],
            vec![1, 1, 1]],
        vec![
            vec![0, 1, 0],
            vec![1, 1, 0]]]
        );

        let mut rotator = Rotator::new(block);
        rotator.rotate_by_i(0);

        assert_eq!(rotator.block.dimensions, Coordinate::new(3, 2, 2));

        rotator.rotate_by_i(1);

        assert_eq!(rotator.block.dimensions, Coordinate::new(2, 2, 3));

        rotator.rotate_by_i(1);
        rotator.rotate_by_i(2);

        assert_eq!(rotator.block.dimensions, Coordinate::new(2, 3, 2));
    }

    #[test]
    fn test_rotation_x_block_body() {
        let block = Block::new(vec![vec![
            vec![0, 1, 1],
            vec![0, 0, 0]],
        vec![
            vec![0, 1, 0],
            vec![1, 1, 0]]]
        );

        let expected_block_after_x_rot = Block::new(vec![vec![
            vec![0, 0, 0],
            vec![1, 1, 0]],vec![
            vec![0, 1, 1],
            vec![0, 1, 0]]]
        );

        let mut rotator = Rotator::new(block);
        rotator.rotate_by_i(0);

        assert_eq!(rotator.block.dimensions, Coordinate::new(3, 2, 2));
        assert_eq!(rotator.block.body.len(), expected_block_after_x_rot.body.len());

        for coord in expected_block_after_x_rot.body {
            assert!(rotator.block.body.contains(&coord));
        }
    }

    #[test]
    fn test_rotation_y_block_body() {
        let block = Block::new(vec![vec![
            vec![0, 1, 1],
            vec![0, 0, 0]],vec![
            vec![0, 1, 0],
            vec![1, 1, 0]]]
        );

        let expected_block_after_y_rot = Block::new(vec![vec![
            vec![1, 0],
            vec![0, 0]],vec![
            vec![1, 1],
            vec![0, 1]],vec![
            vec![0, 0],
            vec![0, 1]]]
        );

        let mut rotator = Rotator::new(block);
        rotator.rotate_by_i(1);

        assert_eq!(rotator.block.dimensions, Coordinate::new(2, 2, 3));
        assert_eq!(rotator.block.body.len(), expected_block_after_y_rot.body.len());

        for coord in expected_block_after_y_rot.body {
            assert!(rotator.block.body.contains(&coord));
        }
    }

    #[test]
    fn test_rotation_z_block_body() {
        let block = Block::new(vec![vec![
            vec![0, 1, 1],
            vec![0, 0, 0]],vec![
            vec![0, 1, 0],
            vec![1, 1, 0]]]
        );

        let expected_block_after_z_rot = Block::new(vec![vec![
            vec![1, 0],
            vec![1, 0],
            vec![0, 0]],vec![
            vec![0, 0],
            vec![1, 1],
            vec![0, 1]]]
        );

        let mut rotator = Rotator::new(block);
        rotator.rotate_by_i(2);

        assert_eq!(rotator.block.dimensions, Coordinate::new(2, 3, 2));
        assert_eq!(rotator.block.body.len(), expected_block_after_z_rot.body.len());

        for coord in expected_block_after_z_rot.body {
            assert!(rotator.block.body.contains(&coord));
        }
    }
}