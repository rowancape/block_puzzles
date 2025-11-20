use std::collections::{HashSet, VecDeque};

use rand::seq::{IndexedMutRandom, IndexedRandom};

use crate::common::*;

pub struct Generator {
    pub field: Field,
    block_unit_sizes: Vec<usize>,
    starting_field_ones: usize
}

#[allow(dead_code)]
impl Generator {
    pub fn new(field: &Field, num_of_blocks: usize, min_block_size: usize) -> Self {
        let starting_field_ones = Self::count_field_ones(&field);
        let free_space = field.volume() - starting_field_ones;
        let block_unit_sizes = Self::divide_into_n_chunks(free_space, num_of_blocks, min_block_size);
        //let block_unit_sizes = vec![free_space / num_of_blocks ; num_of_blocks];

        Self {
            field: field.clone(),
            block_unit_sizes,
            starting_field_ones
        }
    }

    pub fn generate_block(&mut self, block_size: usize, block_index: usize) -> Option<Block> {
        let backup = self.field.clone();
        let mut iterations = 0;
    
        loop {
            self.field = backup.clone();
            let mut field_block = gen_empty_grid(self.field.len_x(), self.field.len_y(), self.field.len_z());
            let mut viable_coords = Vec::new();
            let mut active_coord = self.pick_random_zero_coordinate().expect("Field is full!");
    
            viable_coords.push(active_coord);
            let mut current_size = 1;
    
            while current_size < block_size {
                let coord = match pop_random(&mut viable_coords) {
                    Some(c) => c,
                    None => panic!("No viable coordinates left!"),
                };
    
                self.place_block_cell(&mut field_block, &coord);
                active_coord = coord;
    
                let new_coords = self.get_unvisited_neighbors(active_coord, &viable_coords);
                viable_coords.extend(new_coords);
                current_size += 1;
            }
    
            // Place the final piece
            if let Some(final_coord) = pop_random(&mut viable_coords) {
                self.place_block_cell(&mut field_block, &final_coord);
            } else {
                panic!("No viable coordinates left for final placement!");
            }
    
            if self.is_field_valid_after_placement(block_index) {
                let trimmed = Self::trim_grid(field_block);
                println!("Successfully formed block:");
                print_grid(&trimmed);
                return Some(Block::new(trimmed));
            }

            if iterations >= 40 {
                return None;
            } else {
                iterations += 1;
            }
        }
    }

    pub fn generate_blocks(&mut self) -> Vec<Block> {
        let starting_field = self.field.clone();

        'from_start: loop {
            self.field = starting_field.clone();
            let mut blocks: Vec<Block> = Vec::new();

            for block_size_i in 0..self.block_unit_sizes.len() {
                let size = self.block_unit_sizes[block_size_i];
                let block = self.generate_block(size, block_size_i);
                
                match block {
                    Some(block) => {
                        blocks.push(block);
                    }
                    None => {
                        println!("BLOCK GENERATION GOT STUCK, STARTING OVER!");
                        continue 'from_start;
                    }
                }
            }

            return blocks;
        }
    }

    fn find_first_zero_coordinate(&self) -> Option<Coordinate> {
        for (z, layer) in self.field.iter().enumerate() {
            for (y, row) in layer.iter().enumerate() {
                for (x, point) in row.iter().enumerate() {
                    if *point == 0 {
                        return Some(Coordinate::new(x, y, z))
                    }
                }
            }
        }
        None
    }

    fn pick_random_zero_coordinate(&self) -> Option<Coordinate> {
        let mut rng = rand::rng();
        let mut zero_coords: Vec<Coordinate> = Vec::new();

        for (z, layer) in self.field.iter().enumerate() {
            for (y, row) in layer.iter().enumerate() {
                for (x, point) in row.iter().enumerate() {
                    if *point == 0 {
                        zero_coords.push(Coordinate::new(x, y, z));
                    }
                }
            }
        }
        if let Some(coord) = zero_coords.choose(&mut rng) {
            Some(coord.clone())
        } else {
            None
        }
    }

    fn count_field_ones(field: &Field) -> usize {
        let mut count = 0;
        
        for layer in field.iter() {
            for row in layer {
                for point in row {
                    if *point == 1 {
                        count += 1;
                    }
                }
            }
        }

        count
    }

    pub fn divide_into_n_chunks(num: usize, n: usize, min_size: usize) -> Vec<usize> {
        assert!(num >= n);
        assert!((min_size * n) <= num);

        let mut rng = rand::rng();

        let mut chunks = vec![min_size ; n];
        let difference = num - (min_size * n);

        for _ in 0..difference {
            if let Some(random_chunk) = chunks.choose_mut(&mut rng) {
                *random_chunk += 1;
            }
        }

        chunks
    }

    fn flood_fill(&self, start: Coordinate) -> usize {
        let mut visited = HashSet::new();
        let mut queue = VecDeque::new();
    
        queue.push_back(start.clone());
        visited.insert(start);
    
        while let Some(coord) = queue.pop_front() {
            for axis in 0..3 {
                for dir in [-1, 1] {
                    let mut neighbor = coord.clone();
                    let neighbor_i = coord[axis] as isize + dir;
    
                    if neighbor_i >= 0 && neighbor_i < self.field.len_by_i(axis) as isize {
                        neighbor[axis] = neighbor_i as usize;
    
                        if self.field[neighbor.z][neighbor.y][neighbor.x] == 0 && !visited.contains(&neighbor) {
                            visited.insert(neighbor.clone());
                            queue.push_back(neighbor);
                        }
                    }
                }
            }
        }
    
        visited.len()
    }

    fn place_block_cell(&mut self, grid: &mut Grid, coord: &Coordinate) {
        grid[coord.z][coord.y][coord.x] = 1;
        self.field[coord.z][coord.y][coord.x] = 1;
    }
    
    fn get_unvisited_neighbors(&self, coord: Coordinate, visited: &Vec<Coordinate>) -> Vec<Coordinate> {
        let mut neighbors = Vec::new();
        for axis in 0..3 {
            for dir in [1, -1] {
                let neighbor_i = coord[axis] as isize + dir;
                if neighbor_i < 0 || neighbor_i >= self.field.len_by_i(axis) as isize {
                    continue;
                }
    
                let neighbor_u = neighbor_i as usize;
                if self.field.get_special_by_i(coord.x, coord.y, coord.z, neighbor_u, axis) == 1 {
                    continue;
                }
    
                let mut new_coord = coord.clone();
                new_coord[axis] = neighbor_u;
    
                if !visited.contains(&new_coord) {
                    neighbors.push(new_coord);
                }
            }
        }
        neighbors
    }
    
    fn is_field_valid_after_placement(&self, block_index: usize) -> bool {
        if let Some(coord) = self.find_first_zero_coordinate() {
            let free_space = self.flood_fill(coord);
            let mut remaining_sum = 0;

            for i in 0..block_index + 1 {
                remaining_sum += self.block_unit_sizes[i];
            }

            free_space == self.field.volume() - self.starting_field_ones - remaining_sum
        } else {
            true
        }
    }

    fn trim_grid(grid: Grid) -> Grid {
        let mut min_x = usize::MAX;
        let mut max_x = 0;
        let mut min_y = usize::MAX;
        let mut max_y = 0;
        let mut min_z = usize::MAX;
        let mut max_z = 0;
    
        for z in 0..grid.len() {
            for y in 0..grid[z].len() {
                for x in 0..grid[z][y].len() {
                    if grid[z][y][x] == 1 {
                        min_x = min_x.min(x);
                        max_x = max_x.max(x);
                        min_y = min_y.min(y);
                        max_y = max_y.max(y);
                        min_z = min_z.min(z);
                        max_z = max_z.max(z);
                    }
                }
            }
        }
    
        let mut trimmed = Vec::new();
        for z in min_z..=max_z {
            let mut layer = Vec::new();
            for y in min_y..=max_y {
                let row = grid[z][y][min_x..=max_x].to_vec();
                layer.push(row);
            }
            trimmed.push(layer);
        }
    
        trimmed
    }
}