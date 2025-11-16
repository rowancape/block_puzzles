use std::ops::{Deref, DerefMut};

use super::types::*;
use crate::Coordinate;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Block {
    pub body: Vec<Coordinate>,
    pub dimensions: Coordinate
}

impl Block {
    pub fn new(grid: Grid) -> Self {
        assert!(!grid.is_empty());
        assert!(!grid[0].is_empty());
        assert!(!grid[0][0].is_empty());

        let x = grid.iter().flatten().flatten().collect::<Vec<&u8>>();
        assert!(x.contains(&&1));

        let layer_len = grid[0].len();
        let row_len = grid[0][0].len();

        for layer in grid.iter() {
            if layer.len() != layer_len {
                panic!("Grid layers aren't of uniform size. Grid must be a cuboidal 3D vec.")
            }
            for row in layer {
                if row.len() != row_len {
                    panic!("Grid rows aren't of uniform size. Grid must be a cuboidal 3D vec.")
                }
            }
        }

        let mut coordinates: Vec<Coordinate> = Vec::new();

        for (z, layer) in grid.iter().enumerate() {
            for (y, row) in layer.iter().enumerate() {
                for (x, point) in row.iter().enumerate() {
                    if *point == 1 {
                        coordinates.push(Coordinate::new(x, y, z));
                    }
                }
            }
        }

        let mut x_max = 0;
        let mut y_max = 0;
        let mut z_max = 0;

        for coord in &coordinates {
            if coord.x > x_max { x_max = coord.x }
            if coord.y > y_max { y_max = coord.y }
            if coord.z > z_max { z_max = coord.z }
        }

        Self { 
            body: coordinates,
            dimensions: Coordinate::new(x_max, y_max, z_max)
        }
    }

    pub fn print(&self) {
        println!("{self:?}")
    }
}

impl Deref for Block {
    type Target = Vec<Coordinate>;
    fn deref(&self) -> &Self::Target {
        &self.body
    }
}

impl DerefMut for Block {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.body
    }
}