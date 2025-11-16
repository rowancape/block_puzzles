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
        let dimensions = self.block.dimensions;

        match i {
            0 => { self.rotate_x(dimensions); }
            1 => { self.rotate_y(dimensions); }
            2 => { self.rotate_z(dimensions); }
            _ => { panic!("Invalid rotation index in Rotator function 'rotate_by_i()'!") }
        }
    }

    fn rotate_x(&mut self, max: Coordinate) {
        for coord in &mut self.block.body {
            let y = coord.y;

            coord.y = coord.z;
            coord.z = max.y - y;
        }

        if self.block.dimensions.y != self.block.dimensions.z {
            let temp = self.block.dimensions.y;
            self.block.dimensions.y = self.block.dimensions.z;
            self.block.dimensions.z = temp; 
        }
    }

    fn rotate_y(&mut self, max: Coordinate) {
        for coord in &mut self.block.body {
            let x = coord.x;

            coord.x = coord.z;
            coord.z = max.x - x;
        }

        if self.block.dimensions.x != self.block.dimensions.z {
            let temp = self.block.dimensions.x;
            self.block.dimensions.x = self.block.dimensions.z;
            self.block.dimensions.z = temp; 
        }
    }

    fn rotate_z(&mut self, max: Coordinate) {
        for coord in &mut self.block.body {
            let x = coord.x;

            coord.x = coord.y;
            coord.y = max.x - x;
        }

        if self.block.dimensions.y != self.block.dimensions.x {
            let temp = self.block.dimensions.y;
            self.block.dimensions.y = self.block.dimensions.x;
            self.block.dimensions.x = temp; 
        }
    }
}