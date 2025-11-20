use crate::common::types::{BlockCombination, Block, Coordinate};
use crate::common::rotator::Rotator;
use crate::common::Field;
use itertools::Itertools;

pub struct Solver {
    pub field: Field,
    pub rotators: Vec<Rotator>,
    working_rot_combos: Vec<BlockCombination>,
    pub rotational_combinations_tested: u64,
    pub start_coord_combos_tested: u64
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
            rotational_combinations_tested: 1,
            start_coord_combos_tested: 0
        }
    }

    /// #### Puts all of the internal private methods to use in order to loop through every possible rotational
    /// ####  and positional combination of blocks for the entire puzzle until either:
    /// 
    /// 1. A solution is found and Some() is returned with a tuple containing the blocks in the rotational
    ///     state for which a solution was found, and the starting positions for each of those blocks that
    ///     resulted in the successful solve.
    /// 
    /// 2. The solver has attempted every possible combination of rotational and positional states possible for the
    ///     puzzle and no solution has been found, meaning no solution for the given parameters is possible
    ///     and thus None is returned.
    /// 
    /// ### Note:
    /// 
    /// The time complexity for this is horrendous, as there is an exponential step when computing every possible
    ///   rotational combination for each block, and this is relative to the number of blocks.
    /// 
    /// The solver can check every possibility in a fairly short duration, even with a large field size, so long as
    ///   the number of blocks is <= 4.
    /// At five or more blocks getting to return None may take hours or days assuming none of the blocks are 
    ///   rotationally homogenous (like a perfect cube).
    /// 
    /// This can be improved significantly, but no polynomial-time algorithm for a general solution is currently known,
    ///   since the problem is NP-complete.
    /// 
    /// The fastest method would be using a method called "dancing links" and an algorithm called "algorithm x"
    ///   together simply known as DLX
    pub fn solve(&mut self) -> Option<(Vec<Block>, Vec<Coordinate>)> {
        while self.rotators.last().unwrap().axis_rot_state != [3, 3, 3]  {
    
            for block_combo in &self.working_rot_combos {
    
                let start_coords_each_block = self.find_valid_start_coords(block_combo);
    
                let start_combos_iterator = start_coords_each_block.into_iter().multi_cartesian_product();
    
                for start_combo in start_combos_iterator {
                    self.start_coord_combos_tested += 1;
    
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

    /// Finds all the start coordinates at which any given block alone can be validly placed on the field
    /// 
    /// Then coalesces those valid start coordinates for each block into a Vec<Vec<Coordinate>> where each
    ///   internal Vec holds the valid starting coordinates for the block of the same index in a block combination
    /// 
    /// For example the block at block_combination[1] has it's valid start coordinates placed in
    ///   this_functions_return_vec[1] 
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

    /// Finds starting coordinates for a block which cause the entire block to fit into the dimensions of field.
    /// 
    /// Does NOT assess whether block overlaps with fields internal 1s. Simply checks which start coords
    ///   would be valid for a single block assuming a fully empty field (all points set to 0).
    fn find_start_coord_candidates(&self, block: &Block) -> Vec<Coordinate> {
        let mut single_block_coords: Vec<Coordinate> = Vec::new();

        let field_dimensions = self.field.dimensions();
        let block_dimensions = &block.dimensions;
            
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

    /// Takes a Vec<Coordinate> of start coordinate candidates. These are coordinates from which the block
    ///   does not extend outside field boundaries, but which may still cause the block to overlap with the fields
    ///   internal void structures (a point set to 1 on the fields 3D vec)
    /// 
    /// This simply filters out all candidate coords that cause a block to overlap with the fields internal 1s.
    fn validate_candidate_coords(&self, candidates: Vec<Coordinate>, block: &Block) -> Vec<Coordinate> {
        let mut valid_coords: Vec<Coordinate> = Vec::new();

        for candidate in candidates {
            if !self.does_candidate_overlap_field(&candidate, block) {
                valid_coords.push(candidate.clone());
            }
        }

        valid_coords
    }

    /// Loops over a single blocks coordinates and checks whether placing them on the field
    ///   at a certain candidate starting coordinate would cause overlap with the field itself.
    /// 
    /// Returns true if overlap happens, and false if not.
    fn does_candidate_overlap_field(&self, candidate: &Coordinate, block: &Block) -> bool {
        for coord in block.iter() {
            if self.field[candidate.z + coord.z][candidate.y + coord.y][candidate.x + coord.x] >= 1 {
                return true
            }
        }
        false
    }

    /// Places blocks in the field according to given starting positions.
    /// 
    /// If an overlap between two blocks is observed => returns false.
    /// 
    /// Otherwise the solution was a success and no blocks overlapped => returns true.
    fn is_solve_success(&self, field: &mut Field, blocks: &Vec<Block>, sps: &Vec<Coordinate>) -> bool {
        // sp = start position — a coordinate representing the starting placement position of a block
        // sp[i] corresponds to blocks[i]
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

    /// Finds the next unique rotation of a block, builds the new combinations based on the new unique rotation
    ///   and saves them in working_rot_combos
    /// 
    /// If no unique rotation is found in any of the possible rotations left, the function returns
    ///   and the solve function will notice that all possible rotations of every block have been attempted
    fn compute_new_rotation(&mut self) {
        // If every block has been rotated every possible way -> return.
        // Otherwise save 
        let Some(b) = self.which_block_to_rotate_next() else { return };

        loop {
            self.rotate_block(&b);

            if self.is_current_rotation_unique(&b) {
                
                let new_rotational_combinations: Vec<BlockCombination> = 
                    self.get_combined_blocks(b).into_iter().multi_cartesian_product().collect();
                self.working_rot_combos = new_rotational_combinations;
                
                self.rotational_combinations_tested += self.working_rot_combos.len() as u64;
                return;
            }

            if self.rotators[b].axis_rot_state == [3, 3, 3] {
                break;
            }
        }

        self.compute_new_rotation();
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

    /// Rotates one block based on index without asking questions.
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

    /// Checks and returns whether a block of certain index currently has a unique rotation.
    // This function could really be replaced with a method on the rotator struct if we get nitpicky.
    fn is_current_rotation_unique(&self, block_index: &usize) -> bool {
        !self.rotators[*block_index].previous_rotations.contains(&self.rotators[*block_index].block)
    }

    /// Produces a Vec that contains N number of Vecs where N is the number of blocks in solver.
    /// 
    /// Inside the inner Vecs N-1 of them contain the previous rotational states of a block.
    /// 
    /// And the final one contains only one block, which is the most recently rotated blocks current state.
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