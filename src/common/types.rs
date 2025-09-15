pub use crate::common::{Block, Coordinate};

pub type Grid = Vec<Vec<Vec<u8>>>;
pub type BlockCombination = Vec<Block>;

pub enum Color {
    Red,
    Yellow,
    Green,
    Blue,
}