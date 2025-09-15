use std::ops::{Deref, DerefMut};
use std::fmt;

use super::types::*;

#[derive(Clone, PartialEq, Eq)]
pub struct Field(pub Grid);

impl Field {
    pub fn new(grid: Grid) -> Self {
        assert!(!grid.is_empty());
        assert!(!grid[0].is_empty());
        assert!(!grid[0][0].is_empty());

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

        Self(grid)
    }

    pub fn len_x(&self) -> usize {
        self[0][0].len()
    }

    pub fn len_y(&self) -> usize {
        self[0].len()
    }

    pub fn len_z(&self) -> usize {
        self.len()
    }

    pub fn len_by_i(&self, i: usize) -> usize {
        match i {
            0 => { self.len_x() },
            1 => { self.len_y() },
            2 => { self.len_z() },
            _ => { panic!("Invalid index!") }
        }
    }

    pub fn dimensions(&self) -> Coordinate {
        Coordinate { x: self.len_x(), y: self.len_y(), z: self.len_z() }
    }

    pub fn volume(&self) -> usize {
        self.len_x() * self.len_y() * self.len_z()
    }

    pub fn set_special_by_i(&mut self, x: usize, y: usize, z: usize, special: usize, axis: usize, val: u8) {
        match axis {
            0 => { self[z][y][special] = val },
            1 => { self[z][special][x] = val },
            2 => { self[special][y][x] = val },
            _ => { panic!("Invalid index!"); },
        }
    }

    pub fn get_special_by_i(&self, x: usize, y: usize, z: usize, special: usize, axis: usize) -> u8 {
        match axis {
            0 => { self[z][y][special] },
            1 => { self[z][special][x] },
            2 => { self[special][y][x] },
            _ => { panic!("Invalid index!"); },
        }
    }

    pub fn print(&self) {
        println!("{self:?}")
    }
}

impl fmt::Debug for Field {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "\n")?;
        for (i, layer) in self.iter().enumerate() {
            if i > 0 {
                write!(f, "——————————\n")?;
            }
            for row in layer {
                write!(f, "{:?}\n", row)?;
            }
        }
        Ok(())
    }
}

impl Deref for Field {
    type Target = Grid;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Field {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}