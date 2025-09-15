mod common;
mod generator;
mod solver;

use common::*;
use generator::*;
use solver::*;
use std::time::Instant;

fn main() {
    let start_time = Instant::now();

    let field = Field::new(
        vec![vec![vec![0, 0, 0, 0],  
                    vec![0, 1, 0, 1],
                    vec![0, 0, 0, 0],
                    vec![0, 0, 1, 1]]; 3]);

    let field2 = Field::new(vec![vec![vec![0; 3]; 3]; 3]);


    let mut generator = Generator::new(&field, 3, 4);
    let blocks = generator.generate_blocks();


    let mut solver = Solver::new(field, blocks);

    match solver.solve() {
        Some(solution) => {
            let blocks: BlockCombination = solution.0;
            let coords: Vec<Coordinate> = solution.1;

            let mut field: Vec<Vec<Vec<&str>>> = solver.field
                .iter()
                .map(|v| 
                    v.iter()
                     .map(|inner| inner.iter().map(|_| "⬛").collect())
                     .collect()
                )
                .collect();

            let mut prev_colors: Vec<&str> = Vec::new();
            let colors = [Color::Blue, Color::Green, Color::Red, Color::Yellow];
            let mut color_i = 0;

            for i in 0..blocks.len() {
                let block_color = &colors[color_i % 4];
                color_i += 1;
                produce_solution_field(&mut field, &blocks[i], &coords[i], block_color, &mut prev_colors);
            }

            print_solution(field);

            println!();
            println!();
            println!("Took {:?} to compute.", start_time.elapsed())
        }
        None => {
            println!("Something has gone awry!")
        }
    }
}
