use std::ops::Rem;

pub struct FluidVolume {
	grid_cell_count: (usize, usize),
	// Note that we split the X and Y velocities because we want to use a single solver for f32 rather than having multiple solvers for f32 and (f32, f32).
	x_velocity: Vec<f32>,
	y_velocity: Vec<f32>,
	density: Vec<f32>,

	// These are needed to compute the next steps.  Rather than ask the user to double buffer and worry about memory, we'll take the hit and manage that internally.
	next_density: Vec<f32>,
	divergence: Vec<f32>, // Although this isn't used directly, we pre-allocate a divergence field so we don't have to re-allocate memory every single step.
	projection: Vec<f32>, // This is the divergence free field.  Again, pre-allocated to prevent memory fragmentation.
	next_projection: Vec<f32>, // This is the divergence free field.  Again, pre-allocated to prevent memory fragmentation.
}

impl FluidVolume {
	fn to_idx(&self, x:usize, y:usize) -> usize {
		// Wrap:
		//x.rem_euclid(self.grid_cell_count.0) + (y.rem_euclid(self.grid_cell_count.1))*self.grid_cell_count.0
		// Cap:
		x.max(0).min(self.grid_cell_count.0) + y.max(0).min(self.grid_cell_count.1)*self.grid_cell_count.0
	}

	/// Get center, left, right, top, bottom indices from the given X and Y coordinates.
	/// Respects wrapping rules and boundary conditions.
	fn get_clrtb_idx(&self, x:usize, y:usize) -> (usize, usize, usize, usize, usize) {
		let center_idx = self.to_idx(x, y);
		let left_idx = self.to_idx(x-1, y);
		let right_idx = self.to_idx(x+1, y);
		let top_idx = self.to_idx(x, y-1);
		let bottom_idx = self.to_idx(x, y+1);
		(center_idx, left_idx, right_idx, top_idx, bottom_idx)
	}

	pub fn new(width: usize, height: usize) -> Self {
		FluidVolume {
			grid_cell_count: (width, height),
			//x_velocity_next: (0..(width*height)).map(|_| { (0f32, 0f32) }).collect::<Vec<(f32, f32)>>(),
			x_velocity: (0..(width*height)).map(|_| { 0f32 }).collect::<Vec<f32>>(),
			y_velocity: (0..(width*height)).map(|_| { 0f32 }).collect(),
			density: (0..(width*height)).map(|_| { 0f32 }).collect(),

			next_density: (0..(width*height)).map(|_| { 0f32 }).collect(),
			divergence: (0..(width*height)).map(|_| { 0f32 }).collect(),
			projection: (0..(width*height)).map(|_| { 0f32 }).collect(),
			next_projection: (0..(width*height)).map(|_| { 0f32 }).collect(),
		}
	}

	//
	// Internal accessors that are of no use outside these classes:
	//

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

	// Private accessors based on the above.

	fn get_divergence(&self, x:usize, y:usize) -> f32 {
		self.get_property(x, y, &self.divergence)
	}

	fn set_divergence(&mut self, x:usize, y:usize, value:f32) {
		let idx = self.to_idx(x, y);
		self.divergence[idx] = value;
	}

	fn get_projection(&self, x:usize, y:usize) -> f32 {
		self.get_property(x, y, &self.projection)
	}

	fn set_projection(&mut self, x:usize, y:usize, value:f32) {
		let idx = self.to_idx(x, y);
		self.projection[idx] = value;
	}

	//
	// Public methods:
	//

	// Density:

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

	// Velocity:

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

	pub fn get_velocity(&self, x:usize, y:usize) -> (f32, f32) {
		(self.get_velocity_x(x, y), self.get_velocity_y(x, y))
	}

	/// Interpolate the velocity based on the nearest four points.
	pub fn sample_velocity_x(&self, x: f32, y: f32) -> f32 {
		self.sample_property(x, y, &self.x_velocity)
	}

	pub fn sample_velocity_y(&self, x: f32, y: f32) -> f32 {
		self.sample_property(x, y, &self.y_velocity)
	}

	pub fn sample_velocity(&self, x: f32, y: f32) -> (f32, f32) {
		(self.sample_velocity_x(x, y),self.sample_velocity_y(x, y))
	}

	// Simulation steps: "accumulate forces, diffuse, update constraints" like verlet.

	pub fn step(&mut self, diffusion_rate:f32, delta_time:f32, iterations:usize) {
		self.diffuse(diffusion_rate, iterations);
		self.project(diffusion_rate, iterations);
		self.advect(delta_time);
		self.project(diffusion_rate, iterations);
	}

	/// Perform a step of the fluid density diffusion and write the output to the output set.
	/// We have output as a separate mutable reference so we can double-buffer and switch between
	/// two pre-allocated fluids.
	/// If iterations is zero, this will use the fast solver.
	/// If iterations is more than zero, this will use Gauss-Seidel.
	fn diffuse(&mut self, diffusion_rate:f32, iterations:usize) {
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

					let idx = self.to_idx(x, y);
					self.next_density[idx] = d_current + diffusion_rate*(s_current - d_current);
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
						let (idx, left, right, up, down) = self.get_clrtb_idx(x, y);
						self.next_density[idx] = (self.density[idx] + diffusion_rate * ((self.density[left] + self.density[right] + self.density[up] + self.density[down])*0.25f32)) / (1.0f32 + diffusion_rate);
					}
				}
				self.density[..].copy_from_slice(&self.next_density.as_slice());
			}
		}

		// Finally, copy the next density to this one.
		self.density[..].copy_from_slice(&self.next_density.as_slice());
	}

	/// Use the velocity of a given cell to 'flow' the density from the appropriate cells into the
	/// centers of their new grid locations.  This does _not_ modify velocity.
	/// The velocity is updated by the 'project' step.
	fn advect(&mut self, delta_time:f32) {
		// Note that this updates both next and current density.
		for y in 1..self.grid_cell_count.1 - 1 {
			for x in 1..self.grid_cell_count.0 - 1 {
				// dx/dy are used to select the point on the grid from which our velocity drives.
				let dx = self.get_velocity_x(x, y)*delta_time;
				let dy = self.get_velocity_y(x, y)*delta_time;
				// The outer flow from this should match the inner flow, conservation of mass.
				// So we can, for simplicity, just sample from this x/y - dx/dy.
				let new_density = self.sample_density(x as f32 - dx, y as f32 - dy);
				let idx = self.to_idx(x, y);
				self.next_density[idx] = new_density;
			}
		}

		// Finally, copy the next density to this one.
		self.density[..].copy_from_slice(&self.next_density.as_slice());
	}

	/// Our fluid sim thus far does not conserve mass.
	/// We use another iterative solver to break the system into a curl and a divergence component.
	/// Helmholtz Decomposition: Any vector field can be expressed as the sum of a field that is free of curl and one that is free of divergence.
	/// Can't compute divergence-free directly, so we compute the curl free part and subtract that from the original.
	/// When the steps are complete, we discard the divergence component to maintain mass.
	fn project(&mut self, diffusion_rate:f32, iterations:usize) {
		// First, we compute the divergence of _this_ field based on the velocity,
		// then we perform helmholtz decomposition using the same Gauss-Seidel trick above,
		// then we subtract off the divergence.
		// We use 'output' as a target but we do compute the local divergence on each step.

		let divergence_norm_factor = 1.0f32 / (self.divergence.len() as f32);
		let x_velocity_scale_factor = self.x_velocity.len() as f32;
		let y_velocity_scale_factor = self.y_velocity.len() as f32;

		// Compute divergence for the result field and reset the projection:
		for y in 1..self.grid_cell_count.1 - 1 {
			for x in 1..self.grid_cell_count.0 - 1 {
				let (idx, left, right, up, down) = self.get_clrtb_idx(x, y);
				let h_flow = self.x_velocity[right] - self.x_velocity[left];
				let v_flow = self.y_velocity[down] - self.y_velocity[up];
				self.divergence[idx] = -0.5*(h_flow+v_flow)*divergence_norm_factor;
				self.projection[idx] = 0.0f32;
				//self.next_projection[idx] = 0.0f32; // Not required because we're setting this in the next step.
			}
		}

		// Iteratively solve for the projection field in output.
		for _ in 0..iterations {
			for y in 1..self.grid_cell_count.1 - 1 {
				for x in 1..self.grid_cell_count.0 - 1 {
					let (idx, left, right, up, down) = self.get_clrtb_idx(x, y);
					// Original paper by Stam has this as the mean of the projection AND divergence, but I think that's wrong, so I'm using the sum of divergence and the mean of projections.
					self.next_projection[idx] = (self.divergence[idx] + diffusion_rate * ((self.projection[left] + self.projection[right] + self.projection[up] + self.projection[down]) * 0.25f32)) / (1.0f32 + diffusion_rate);
				}
			}
			self.projection[..].copy_from_slice(&self.next_projection.as_slice());
		}

		// Finally update our velocity based on the iterated projection.
		for y in 1..self.grid_cell_count.1 - 1 {
			for x in 1..self.grid_cell_count.0 - 1 {
				let (idx, left, right, up, down) = self.get_clrtb_idx(x, y);
				self.x_velocity[idx] -= 0.5*(self.projection[right]-self.projection[left])*x_velocity_scale_factor;
				self.y_velocity[idx] -= 0.5*(self.projection[down]-self.projection[up])*y_velocity_scale_factor;
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
	fn test_idx_access() {
		let fluid = FluidVolume::new(3, 3);
		// 0, 1, 2
		// 3, 4, 5
		// 6, 7, 8
		let (center, left, right, up, down) = fluid.get_clrtb_idx(1, 1);
		assert_eq!(center, 4);
		assert_eq!(left, 3);
		assert_eq!(right, 5);
		assert_eq!(up, 1);
		assert_eq!(down, 7);
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
}