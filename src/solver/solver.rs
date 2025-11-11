use crate::common::types::{BlockCombination, Block, Coordinate};
use crate::common::rotator::Rotator;
use crate::common::Field;
use itertools::Itertools;

pub struct Solver {
    pub field: Field,
    pub rotators: Vec<Rotator>,
    working_rot_combos: Vec<BlockCombination>
}

impl Solver {
    pub fn new(field: Field, blocks: Vec<Block>) -> Solver {
        let working_rot_combos: Vec<BlockCombination> = 
            vec![blocks.iter().map(|block| block.clone()).collect()];

        let mut rotators: Vec<Rotator> = Vec::new();
        for block in blocks {
            rotators.push(Rotator::new(block));
        }

        Solver {
            field,
            rotators,
            working_rot_combos,
        }
    }

    pub fn solve(&mut self) -> Option<(Vec<Block>, Vec<Coordinate>)> {
        while self.rotators.last().unwrap().axis_rot_state != [3, 3, 3]  {
    
            for block_combo in &self.working_rot_combos {
    
                let start_coords_each_block = self.find_valid_start_coords(block_combo);
    
                let start_combos_iterator = start_coords_each_block.into_iter().multi_cartesian_product();
    
                for start_combo in start_combos_iterator {
    
                    let mut active_field = self.field.clone();
    
                    if self.is_solve_success(&mut active_field, block_combo, &start_combo) {
                        return Some((block_combo.clone(), start_combo));
                    }
                }
            }
            self.compute_new_rotation();
        }

        None
    }

    fn find_valid_start_coords(&self, combo: &BlockCombination) -> Vec<Vec<Coordinate>> {
        // Possible start coordinates for each block in a combination of blocks.
        // If there are four blocks in a combination this would be vector of four vectors,
        // each containing some number of possible start coordinates for the corresponding block.
        let mut start_coords_each_block: Vec<Vec<Coordinate>> = Vec::new();

        for block in combo {
            let candidates = self.find_start_coord_candidates(block);
            let valid_start_coords = self.validate_candidate_coords(candidates, block);

            start_coords_each_block.push(valid_start_coords);
        }

        start_coords_each_block
    }

    fn find_start_coord_candidates(&self, block: &Block) -> Vec<Coordinate> {
        let mut single_block_coords: Vec<Coordinate> = Vec::new();

        let field_dimensions = self.field.dimensions();
        let block_dimensions = block.dimensions();
            
        for z in 0..field_dimensions.z {
            // If block can't fit in the z axis, break.
            if field_dimensions.z - z < block_dimensions.z + 1 {
                break;
            }

            for y in 0..field_dimensions.y {
                // If block can't fit in the y axis, break.
                if field_dimensions.y - y < block_dimensions.y + 1 {
                    break;
                }

                for x in 0..field_dimensions.x {
                    // If block can't fit in the x axis, break.
                    if field_dimensions.x - x < block_dimensions.x + 1 {
                        break;
                    }
                    // If block fits at start coordinate [x, y, z] in all axis then save that coord.
                    single_block_coords.push(Coordinate::new(x, y, z));
                }
            }
        }

        single_block_coords
    }

    fn validate_candidate_coords(&self, candidates: Vec<Coordinate>, block: &Block) -> Vec<Coordinate> {
        let mut valid_coords: Vec<Coordinate> = Vec::new();

        for candidate in candidates {
            if !self.does_candidate_overlap_field(candidate, block) {
                valid_coords.push(candidate);
            }
        }

        valid_coords
    }

    fn does_candidate_overlap_field(&self, candidate: Coordinate, block: &Block) -> bool {
        for coord in block.iter() {
            if self.field[candidate.z + coord.z][candidate.y + coord.y][candidate.x + coord.x] >= 1 {
                return true
            }
        }
        false
    }

    fn is_solve_success(&self, field: &mut Field, blocks: &Vec<Block>, sps: &Vec<Coordinate>) -> bool {
        for (spi, sp) in sps.iter().enumerate() {
            for block_coord in blocks[spi].iter() {
                if field[block_coord.z + sp.z][block_coord.y + sp.y][block_coord.x + sp.x] >= 1 {
                    return false;
                }
                else {
                    field[block_coord.z + sp.z][block_coord.y + sp.y][block_coord.x + sp.x] += 1
                }
            }
        }
        
        true
    }

    fn compute_new_rotation(&mut self) {
        if let Some(b) = self.which_block_to_rotate_next() {
            
            loop {
                self.rotate_block(&b);
                if self.rotators[b].axis_rot_state == [3, 3, 3] {
                    break;
                }
                else if self.is_current_rotation_unique(&b) {
                    let unique_block = self.rotators[b].block.clone();
                    self.rotators[b].previous_rotations.insert(unique_block);
                    break;
                }
            }

            let temp = self.get_combined_blocks(b);
            let new_rotational_combinations: Vec<BlockCombination> = temp.into_iter().multi_cartesian_product().collect();
            self.working_rot_combos = new_rotational_combinations;
        }
    }

    /// Loops over the rotators and returns either:
    /// - Some(i) with the index of the rotator that should be rotated next.
    /// - None if every rotation of all rotators have been tried already (no more new rotations to try).
    /// 
    /// If a rotators axis_rot_state is at [3, 3, 3], that means every possible rotation has been tried already.
    fn which_block_to_rotate_next(&self) -> Option<usize> {
        for (i, rotator) in self.rotators.iter().enumerate() {
            if rotator.axis_rot_state == [3, 3, 3] {
                continue;
            }
            else {
                return Some(i);
            }
        } 

        None
    }


    fn rotate_block(&mut self, block_index: &usize) {
        let rotator = &mut self.rotators[*block_index];

        for i in 0..3 {
            if rotator.axis_rot_state[i] == 3 {
                rotator.rotate_by_i(i);
                continue;
            }
            else {
                rotator.rotate_by_i(i);
                break;
            }
        }
    }

    fn is_current_rotation_unique(&self, block_index: &usize) -> bool {
        !self.rotators[*block_index].previous_rotations.contains(&self.rotators[*block_index].block)
    }

    fn get_combined_blocks(&self, b: usize) -> Vec<BlockCombination> {
        let mut result: Vec<BlockCombination> = Vec::new();

        for (i, rotator) in self.rotators.iter().enumerate() {
            if i == b {
                result.push(vec![rotator.block.clone()]);
            } 
            else {
                result.push(rotator.previous_rotations.iter().cloned().collect());
            }
        }

        result
    }
}