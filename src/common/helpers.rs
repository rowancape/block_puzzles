use rand::Rng;

use super::types::*;

pub fn print_grid(grid: &Grid) {
    println!();
    for (i, layer) in grid.iter().enumerate() {
        if i != 0 {
            println!("——————————");
        }
        for row in layer {
            println!("{row:?}");
        }
    }
    println!();
}

pub fn print_grid_from_coords(coords: &Vec<Coordinate>) {
    let max_x = coords.iter().map(|c| c.x).max().unwrap();
    let max_y = coords.iter().map(|c| c.y).max().unwrap();
    let max_z = coords.iter().map(|c| c.z).max().unwrap();

    let mut grid = vec![vec![vec![0u8; max_x + 1]; max_y + 1]; max_z + 1];

    for coord in coords {
        grid[coord.z][coord.y][coord.x] = 1;
    }

    print_grid(&grid);
}

pub fn pop_random<T>(vec: &mut Vec<T>) -> Option<T> {
    if vec.is_empty() {
        None
    } else {
        let mut rng = rand::rng();
        let index = rng.random_range(0..vec.len());
        Some(vec.swap_remove(index))
    }
}

pub fn gen_empty_grid(x: usize, y: usize, z: usize) -> Grid {
    vec![vec![vec![0; x]; y]; z]
}

pub fn print_solution(field: Vec<Vec<Vec<&str>>>) {
    for layer in field {
        println!();
        for row in layer {
            println!();
            for point in row {
                print!("{point}");
            }
        }
    }
}

pub fn produce_solution_field(
    field: &mut Vec<Vec<Vec<&str>>>,
    block: &Block,
    coord: &Coordinate,
    color: &Color,
    prev_colors: &mut Vec<&str>,
) {
    let mut c = "";

    match color {
        Color::Blue => {
            if prev_colors.contains(&"🟦") {
                c = "🔵"
            } else {
                c = "🟦"
            }
        }
        Color::Red => {
            if prev_colors.contains(&"🟥") {
                c = "🔴"
            } else {
                c = "🟥"
            }
        }
        Color::Yellow => {
            if prev_colors.contains(&"🟨") {
                c = "🟡"
            } else {
                c = "🟨"
            }
        }
        Color::Green => {
            if prev_colors.contains(&"🟩") {
                c = "🟢"
            } else {
                c = "🟩"
            }
        }
    }

    prev_colors.push(c);

    for b_coord in block.iter() {
        field[coord.z + b_coord.z][coord.y + b_coord.y][coord.x + b_coord.x] = c;
    }
}
