use std::ops::Rem;

pub struct FluidVolume {
	grid_cell_count: (usize, usize),
	// Note that we split the X and Y velocities because we want to use a single solver for f32 rather than having multiple solvers for f32 and (f32, f32).
	x_velocity: Vec<f32>,
	y_velocity: Vec<f32>,
	density: Vec<f32>,
}

impl FluidVolume {
	#[inline]
	fn to_idx(&self, x:usize, y:usize) -> usize {
		x.rem_euclid(self.grid_cell_count.0) + (y.rem_euclid(self.grid_cell_count.1))*self.grid_cell_count.0
	}

	pub fn new(width: usize, height: usize) -> Self {
		FluidVolume {
			grid_cell_count: (width, height),
			//x_velocity_next: (0..(width*height)).map(|_| { (0f32, 0f32) }).collect::<Vec<(f32, f32)>>(),
			x_velocity: (0..(width*height)).map(|_| { 0f32 }).collect::<Vec<f32>>(),
			y_velocity: (0..(width*height)).map(|_| { 0f32 }).collect(),
			density: (0..(width*height)).map(|_| { 0f32 }).collect(),
		}
	}

	fn get_property(&self, x:usize, y:usize, field: &Vec<f32>) -> f32 {
		field[self.to_idx(x, y)]
	}

	fn sample_property(&self, x: f32, y: f32, field: &Vec<f32>) -> f32 {
		let (x_cell, x_frac) = get_fractional_cell_component(x, self.grid_cell_count.0);
		let (y_cell, y_frac) = get_fractional_cell_component(y, self.grid_cell_count.1);

		let top = lerp(
			self.get_property(x_cell, y_cell, field),
			self.get_property(x_cell+1, y_cell, field),
			x_frac
		);

		let bottom = lerp(
			self.get_property(x_cell, y_cell+1, field),
			self.get_property(x_cell+1, y_cell+1, field),
			x_frac
		);

		lerp(top, bottom, y_frac)
	}

	/// Horizontally and vertically interpolate a measurement from the surrounding grid squares.
	/// Assumes that x is somewhere between 0 and the number of cells wide.
	/// Assumes y is also somewhere between 0 and the number of cells high.
	pub fn sample_density(&self, x: f32, y: f32) -> f32 {
		self.sample_property(x, y, &self.density)
	}

	/// Return the density of the given cell, no interpolation.
	pub fn get_density(&self, x:usize, y:usize) -> f32 {
		self.get_property(x, y, &self.density)
	}

	pub fn set_density(&mut self, x:usize, y:usize, value:f32) {
		let idx = self.to_idx(x, y);
		self.density[idx] = value;
	}

	pub fn set_velocity(&mut self, x:usize, y:usize, velocity:(f32, f32)) {
		let idx = self.to_idx(x, y);
		self.x_velocity[idx] = velocity.0;
		self.y_velocity[idx] = velocity.1;
	}

	pub fn get_velocity_x(&self, x:usize, y:usize) -> f32 {
		self.get_property(x, y, &self.x_velocity)
	}

	pub fn get_velocity_y(&self, x:usize, y:usize) -> f32 {
		self.get_property(x, y, &self.y_velocity)
	}

	/// Interpolate the velocity based on the nearest four points.
	pub fn sample_velocity_x(&self, x: f32, y: f32) -> f32 {
		self.sample_property(x, y, &self.x_velocity)
	}

	pub fn sample_velocity_y(&self, x: f32, y: f32) -> f32 {
		self.sample_property(x, y, &self.y_velocity)
	}

	/// Perform a step of the fluid density diffusion and write the output to the output set.
	/// We have output as a separate mutable reference so we can double-buffer and switch between
	/// two pre-allocated fluids.
	/// If iterations is zero, this will use the fast solver.
	/// If iterations is more than zero, this will use Gauss-Seidel.
	pub fn step_density_diffusion(&self, k:f32, iterations:usize, output: &mut Self) {
		// d(x,y) = density at x,y
		// s(x,y) = mean of all adjacent squares to x,y = (d(x+1,y) + d(x-1,y) + d(x,y+1) + d(x,y-1))/4.
		// k = rate of change

		// Could use:
		// d_next = d_current + k * (s_current - d_current) // Lerp from current to average of adjacent.
		// But this is not stable.

		if iterations == 0 {
			for y in 1..self.grid_cell_count.1 {
				for x in 1..self.grid_cell_count.0 {
					let d_current = self.get_density(x, y);
					let s_current = (self.get_density(x-1, y) + self.get_density(x+1, y) + self.get_density(x, y-1) + self.get_density(x, y+1))*0.25f32;
					output.set_density(x, y, d_current + k*(s_current - d_current));
				}
			}
		}

		// Reframing and making the equation hyperbolic for better stability:
		// d_current = d_next - k * (s_next - d_next)
		// d_current = d_next - k*s_next + k*d_next
		// d_current - k*s_next = d_next + k*d_next
		// d_current - k*s_next = d_next * (1 + k)
		// (d_current - k*s_next)/(1 + k) = d_next
		// d_next(x,y) = (d_current(x,y) + k*(mean(d_next(...))) / 1+k
		// We can solve with a few iterations of Gauss-Seidel
		if iterations > 0 {
			for _ in 0..iterations {
				for y in 1..self.grid_cell_count.1 - 1 {
					for x in 1..self.grid_cell_count.0 - 1 {
						output.set_density(x, y, (
							self.get_density(x, y) + k*(output.get_density(x-1, y) + output.get_density(x+1, y) + output.get_density(x,y-1) + output.get_density(x, y+1))*0.25f32) /
							(1.0f32 + k)
						);
					}
				}
			}
		}
	}

	pub fn step_velocity(&self, delta_time:f32, output: &mut Self) {
		for y in 1..self.grid_cell_count.1 - 1 {
			for x in 1..self.grid_cell_count.0 - 1 {
				// dx/dy are used to select the point on the grid from which our velocity drives.
				let dx = self.get_velocity_x(x, y);
				let dy = self.get_velocity_y(x, y);
				// Sample from the velocities at this point.
				// We need to cap + wrap these, unlike the other cases.
				output.set_velocity(x, y, (0.0, 0.0));
			}
		}
	}
}


/// Linearly interpolate between a and b by the amount 'amount'.
/// Example:
/// ```rust
/// assert_eq!(lerp(0.0f32, 1.0f32, 0.5), 0.5);
///
/// assert_eq!(lerp(0.0f32, 100.0f32, 0.1), 10.0f32);
/// assert_eq!(lerp(0.0f32, 100.0f32, 0.5), 50.0f32);
/// assert_eq!(lerp(0.0f32, 100.0f32, 0.9), 90.0f32);
///
/// assert_eq!(lerp(50.0f32, 100.0f32, 0.5), 75.0f32);
/// ```
fn lerp(a:f32, b:f32, amount:f32) -> f32 { a + (b-a)*amount }


/// Given a grid size and a floating position, return
/// a) the real part that represents the floor of the cell and
/// b) the fractional component that represents the position between this and the next cell.
/// We assume that each cell in the grid is 1 unit wide for consistency with the get_#_methods.
///
/// Example:
/// ```rust
/// let cell, frac = get_fractional_cell_component(1.4, 4);
/// assert_eq!(cell, 1);
/// assert_eq!(frac, 0.15);
/// ```
fn get_fractional_cell_component(pos: f32, cell_count:usize) -> (usize, f32) {
	let cell_idx:usize = (pos.round() as usize).rem_euclid(cell_count);
	let fraction:f32 = pos.fract();
	(cell_idx, fraction)
}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn test_fractional_cell_component() {
		assert_eq!(get_fractional_cell_component(5.0f32, 10), (5, 0.0));
	}

	#[test]
	fn test_boundary_conditions() {
		let mut fluid = FluidVolume::new(3, 1);
		fluid.set_density(0, 0, 1.0);
		fluid.set_density(1, 0, 2.0);
		fluid.set_density(2, 0, 3.0);
		// Test Wrap:
		assert_eq!(fluid.get_density(3, 0), 1f32);
		assert_eq!(fluid.get_density(4, 0), 2f32);
		assert_eq!(fluid.get_density(5, 0), 3f32);
	}

	#[test]
	fn sanity_test_density() {
		let mut fluid = FluidVolume::new(10, 10);
		let mut next_fluid = FluidVolume::new(10, 10);
		fluid.set_density(5, 5, 1.0);
		fluid.step_density_diffusion(0.1, 0, &mut next_fluid);
		next_fluid.step_density_diffusion(0.1, 0, &mut fluid);
		fluid.step_density_diffusion(0.1, 0, &mut next_fluid);
		next_fluid.step_density_diffusion(0.1, 0, &mut fluid);
		fluid.step_density_diffusion(0.1, 0, &mut next_fluid);
		next_fluid.step_density_diffusion(0.1, 0, &mut fluid);
	}
}