mod common;
mod generator;
mod solver;
mod dlx;

use common::*;
use generator::*;
use solver::*;
use dlx::*;
use std::time::Instant;

fn main() {
    let start_time = Instant::now();

    let field = Field::new(
        vec![vec![vec![0, 0, 0, 1],  
                        vec![1, 0, 0, 0],
                        vec![0, 0, 0, 0]]; 2]);
                    
    let block1 = Block::new(
        vec![vec![
                vec![1, 1],  
                vec![0, 1],
                vec![0, 1]],
            vec![
                vec![0, 0],  
                vec![0, 0],
                vec![0, 1]]]);

    let block2 = Block::new(
        vec![vec![
                vec![1, 1],  
                vec![0, 1],
                vec![1, 1]]]);

    let block3 = Block::new(
        vec![vec![
                vec![0, 1],  
                vec![1, 1],
                vec![0, 1]],
            vec![
                vec![0, 1],  
                vec![0, 0],
                vec![0, 0]]]);

    let block4 = Block::new(
        vec![vec![
                vec![1, 1],  
                vec![1, 1]],
            vec![
                vec![0, 1],  
                vec![0, 0]]]);

    let blocks = vec![block1, block2, block3, block4];

    let mut dlx = DLX::new(field, blocks);
    let solutions = dlx.search();

    println!("TIME: {:?}", start_time.elapsed());

    for solution in solutions {
        println!("\n{:?}", solution);
    }

    // let mut generator = Generator::new(&field, 4, 5);
    // let blocks = generator.generate_blocks();

    /*
    let mut solver = Solver::new(field, blocks);

    match solver.solve() {
        Some((blocks, coords)) => {

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
            println!("Took {:?} to compute.", start_time.elapsed());
            println!();
            println!("Number of block position combinations tested: {}", solver.start_coord_combos_tested);
            println!();
            println!("Number of block rotation state combinations tested: {}", solver.rotational_combinations_tested);
        }
        None => {
            println!("No solution found!");
            println!();
            println!("Took {:?} to go through all possible arrangments.", start_time.elapsed());
            println!();
            println!("Number of block position combinations tested: {}", solver.start_coord_combos_tested);
            println!();
            println!("Number of block rotation state combinations tested: {}", solver.rotational_combinations_tested);
        }
    }
    */
}
