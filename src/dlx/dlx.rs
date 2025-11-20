use crate::common::*;

#[allow(dead_code)]
pub struct DLX {
    // header/root entry node
    root: usize,
    // Total number of nodes in the matrix.
    n_nodes: usize,
    // number of columns and rows in total
    n_rows: usize,
    n_cols: usize,
    // n_nodes_in_col[i] contains the # of nodes in column with index i. 
    n_nodes_in_col: Vec<usize>,
    // row[i]/col[i] contains the row/col that the node of index i is in.
    row: Vec<usize>,
    col: Vec<usize>,
    // left[i] contains the index of the node that is to the left of the node with index i and so on.
    left: Vec<usize>,
    right: Vec<usize>,
    up: Vec<usize>,
    down: Vec<usize>,
}

pub struct DlxRow {
    cells_filled: Vec<usize>,
    block_index: usize
}


#[allow(dead_code)]
impl DLX {
    pub fn new(field: Field, blocks: Vec<Block>) -> DLX {
        let n_blocks = blocks.len();
        let rotators: Vec<Rotator> = blocks.into_iter().map(|block| Rotator::new(block)).collect();
        let mut dlx_rows: Vec<DlxRow> = Vec::new();

        let mut n_rows = 1;
        let n_cols = field.get_n_field_zeros() + n_blocks;
        // Number of nodes is set to n_cols for header row and root node (+ 1)
        let mut n_nodes = n_cols + 1;
        let mut n_nodes_in_col: Vec<usize> = vec![0; n_cols];
        
        // For each block
        for (rot_i, mut rotator) in rotators.into_iter().enumerate() {
            // While the block hasn't been rotated every possible way yet
            while rotator.axis_rot_state != [3, 3, 3] {
                // Rotate to the next potentially unique rotation
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
                // If it's not unique, go to the next iteration of the while loop
                if !rotator.is_current_rot_unique() {
                    continue;
                }

                let start_coords = Self::starting_positions(&field, &rotator.block);
                let mut valid_coords = Vec::new();

                for sc in start_coords {
                    if Self::is_start_coord_valid(&field, &rotator.block.body, &sc) {
                        valid_coords.push(sc);
                    }
                }

                n_rows += valid_coords.len();
                for sc in valid_coords {
                    let filled_cells = Self::which_cells_filled(&field, &rotator.block.body, &sc);
                    n_nodes += filled_cells.len() + 1;
                    dlx_rows.push(DlxRow { cells_filled: filled_cells, block_index: rot_i });
                }
            }
        }

        let mut dlx_cols: Vec<Vec<usize>> = vec![vec![]; n_cols];

        let mut left: Vec<usize> = vec![0; n_nodes];
        let mut right: Vec<usize> = vec![0; n_nodes];
        let mut up: Vec<usize> = vec![0; n_nodes];
        let mut down: Vec<usize> = vec![0; n_nodes];
        let mut row: Vec<usize> = vec![0; n_nodes];
        let mut col: Vec<usize> = vec![0; n_nodes];

        let root: usize = 0;
        left[root] = n_cols;
        right[root] = 1;

        for header_node in 1..=n_cols {
            dlx_cols[header_node - 1].push(header_node);
            col[header_node] = header_node - 1;
            row[header_node] = 0;

            left[header_node] = header_node - 1;
            
            // If last header node
            if header_node == n_cols {
                right[root];
            } else {
                right[header_node] = header_node + 1;
            }
        }

        let mut current_node = n_cols + 1;
        for (row_i, unprocessed_row) in dlx_rows.iter().enumerate() {
            row[current_node] = row_i + 1;
            col[current_node] = unprocessed_row.block_index;
            left[current_node] = unprocessed_row.cells_filled.len() + current_node;
            right[current_node] = current_node + 1;
            dlx_cols[unprocessed_row.block_index].push(current_node);
            
            current_node += 1;
            
            for (i, nth_field_cell) in unprocessed_row.cells_filled.iter().enumerate() {
                let true_col = n_blocks + nth_field_cell;
                row[current_node] = row_i + 1;
                col[current_node] = true_col;
                dlx_cols[true_col].push(current_node);
                left[current_node] = current_node - 1;

                // If last node in row
                if i != unprocessed_row.cells_filled.len() - 1 {
                    right[current_node] = current_node - unprocessed_row.cells_filled.len();
                } else {
                    right[current_node] = current_node + 1;
                }

                current_node += 1;
            }
        }

        for (col_i, col) in dlx_cols.iter().enumerate() {
            for (node_i, node) in col.iter().enumerate() {
                n_nodes_in_col[col_i] += 1;
                // If first node in col
                if node_i == 0 {
                    up[*node] = *col.last().unwrap()
                } else {
                    up[*node] = col[node_i - 1];
                }
                
                // If last node in col
                if node_i == col.len() - 1 {
                    down[*node] = col[0];
                } else {
                    down[*node] = col[node_i + 1];
                }
            }
        }

        Self {
            root,
            n_nodes,
            n_rows,
            n_cols,
            n_nodes_in_col,
            row,
            col,
            left,
            right,
            up,
            down
        }
    }

    pub fn print_data(&self) {
        println!("\nNUMBER OF NODES IN SPARSE MATRIX: {}\n", self.n_nodes);

        for i in 0..self.n_cols {
            println!("NODES IN COL {i}: {}\n", self.n_nodes_in_col[i]);
        }
        println!("NUMBER OF ROWS: {}\n", self.n_rows);
        println!("RIGHT OF ROOT (0): {}", self.right[0]);
        println!("DOWN OF FIRST FIELD CELL COL: {}", self.down[5]);

        println!("WHAT IS BELOW 25: {}", self.down[25]);

        let col_header = 2;
        let mut curr_node_in_col = self.down[col_header];
        
        while col_header != curr_node_in_col {
            println!("NODE: {}, ROW: {}", curr_node_in_col, self.row[curr_node_in_col]);
            curr_node_in_col = self.down[curr_node_in_col];
        }
    }

    fn is_start_coord_valid(field: &Field, block: &Vec<Coordinate>, sc: &Coordinate) -> bool {
        
        for bc in block {
            if field[sc.z+bc.z][sc.y+bc.y][sc.x+bc.x] == 1 {
                return false;
            }
        }

        true
    }

    fn which_cells_filled(field: &Field, block: &Vec<Coordinate>, sc: &Coordinate) -> Vec<usize> {
        let mut cells_filled = Vec::new();
        let mut cell_count: usize = 0;

        let block: Vec<Coordinate> = block.iter()
            .map(|b| Coordinate::new(b.x + sc.x, b.y + sc.y, b.z + sc.z)).collect();

        for (z, layer) in field.iter().enumerate() {
            for (y, row) in layer.iter().enumerate() {
                for (x, point) in row.iter().enumerate() {
                    if block.contains(&Coordinate::new(x, y, z)) {
                        cells_filled.push(cell_count);
                    }
                    if *point == 0 {
                        cell_count += 1;
                    }
                }
            }
        }

        cells_filled
    }

    fn starting_positions(field: &Field, block: &Block) -> Vec<Coordinate> {
        let mut result: Vec<Coordinate> = Vec::new();

        let field_dimensions = field.dimensions();
        let block_dimensions = &block.dimensions;
            
        for z in 0..field_dimensions.z {
            // If block can't fit in the z axis, break.
            if field_dimensions.z - z < block_dimensions.z {
                break;
            }

            for y in 0..field_dimensions.y {
                // If block can't fit in the y axis, break.
                if field_dimensions.y - y < block_dimensions.y {
                    break;
                }

                for x in 0..field_dimensions.x {
                    // If block can't fit in the x axis, break.
                    if field_dimensions.x - x < block_dimensions.x {
                        break;
                    }
                    // If block fits at start coordinate [x, y, z] in all axis then save that coord.
                    result.push(Coordinate::new(x, y, z));
                }
            }
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_which_cells_filled_basic() {
        // 3×3×1 field with some filled cells
        // Layer 0:
        // 0 1 0
        // 0 0 0
        // 1 0 0
        let field = Field::new(vec![vec![
            vec![0, 1, 0],
            vec![0, 0, 0],
            vec![1, 0, 0]]]
        );

        // Block: 2×2 square (relative coords)
        let block = vec![
            Coordinate::new(0, 0, 0),
            Coordinate::new(1, 0, 0),
            Coordinate::new(0, 1, 0),
            Coordinate::new(1, 1, 0),
        ];

        // Place block at start coordinate (1, 1, 0)
        let start = Coordinate::new(1, 1, 0);

        let filled = DLX::which_cells_filled(&field, &block, &start);

        assert_eq!(filled, vec![3, 4, 5, 6]);
    }

    #[test]
    fn test_which_cells_filled_3d() {
        // 2×2×2 field, all empty
        let field = Field::new(vec![vec![
            vec![0, 0],
            vec![0, 0]],vec![
            vec![0, 0],
            vec![0, 0]]]
        );

        let block = Block::new(vec![vec![
            vec![0, 1],
            vec![1, 1]],vec![
            vec![0, 0],
            vec![0, 1]]]
        );

        let start = Coordinate::new(0, 0, 0);

        let filled = DLX::which_cells_filled(&field, &block, &start);

        assert_eq!(filled, vec![1, 2, 3, 7]);
    }

    #[test]
    fn test_is_start_coord_valid() {
        let field = Field::new(vec![vec![
            vec![0, 0, 0],
            vec![0, 0, 0],
            vec![1, 0, 0]],
        vec![
            vec![0, 0, 0],
            vec![0, 0, 1],
            vec![1, 0, 0]]]
        );

        let block = Block::new(vec![vec![
            vec![0, 1],
            vec![1, 1]],
        vec![
            vec![0, 0],
            vec![0, 1]]]
        );

        let start1 = Coordinate::new(0, 0, 0);
        let start2 = Coordinate::new(1, 0, 0);
        let start3 = Coordinate::new(0, 1, 0);
        let start4 = Coordinate::new(1, 1, 0);

        assert!(DLX::is_start_coord_valid(&field, &block, &start1));
        assert!(!DLX::is_start_coord_valid(&field, &block, &start2));
        assert!(!DLX::is_start_coord_valid(&field, &block, &start3));
        assert!(DLX::is_start_coord_valid(&field, &block, &start4));
    }

    #[test]
    #[should_panic]
    fn test_is_start_coord_valid_fail() {
        let field = Field::new(vec![vec![
            vec![0, 0, 0],
            vec![0, 0, 0],
            vec![1, 0, 0]],
        vec![
            vec![0, 0, 0],
            vec![0, 0, 1],
            vec![1, 0, 0]]]
        );

        let block = Block::new(vec![vec![
            vec![0, 1],
            vec![1, 1]],vec![
            vec![0, 0],
            vec![0, 1]]]
        );

        let start1 = Coordinate::new(1, 0, 0);

        assert!(DLX::is_start_coord_valid(&field, &block, &start1));
    }

    #[test]
    fn test_starting_positions() {
        let field = Field::new(vec![vec![
            vec![1, 0],
            vec![0, 0]]; 2]
        );

        let block = Block::new(vec![vec![
            vec![1, 1],
            vec![1, 1]]; 2]
        );

        let start1 = Coordinate::new(0, 0, 0);

        let calculated_coords = DLX::starting_positions(&field, &block);

        assert_eq!(calculated_coords[0], start1);
        assert_eq!(calculated_coords.len(), 1);
    }

    #[test]
    fn test_starting_positions_multiple() {
        let field = Field::new(vec![vec![
            vec![1, 0, 1],
            vec![0, 0, 0],
            vec![0, 0, 0]]; 3]
        );

        let block = Block::new(vec![vec![
            vec![1],
            vec![1]]; 2]
        );

        let expected_coords = vec![
            Coordinate::new(0, 0, 0),
            Coordinate::new(1, 0, 0),
            Coordinate::new(2, 0, 0),
            Coordinate::new(0, 0, 1),
            Coordinate::new(1, 0, 1),
            Coordinate::new(2, 0, 1),
            Coordinate::new(0, 1, 0),
            Coordinate::new(1, 1, 0),
            Coordinate::new(2, 1, 0),
            Coordinate::new(0, 1, 1),
            Coordinate::new(1, 1, 1),
            Coordinate::new(2, 1, 1),
        ];

        let calculated_coords = DLX::starting_positions(&field, &block);

        assert_eq!(expected_coords.len(), calculated_coords.len());

        for coordinate in expected_coords {
            assert!(calculated_coords.contains(&coordinate));
        }
    }
}