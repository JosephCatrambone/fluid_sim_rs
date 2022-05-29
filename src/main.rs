
mod fluid;

use minifb::{Key, MouseButton, MouseMode, Window, WindowOptions};
use crate::fluid::FluidVolume;

const WIDTH: usize = 480;
const HEIGHT: usize = 480;
const DOWNSCALE: usize = 5;
const DIFFUSION_RATE: f32 = 0.5f32;
const SUBSTEPS: usize = 2;
const SIMULATION_ITERATIONS: usize = 0;
const INFLOW_RATE: f32 = 5.0f32;

fn main() {
	let mut buffer: Vec<u32> = vec![0; WIDTH * HEIGHT];
	let mut fluid_tick = fluid::FluidVolume::new(WIDTH/DOWNSCALE, HEIGHT/DOWNSCALE);
	let mut fluid_tock = fluid::FluidVolume::new(WIDTH/DOWNSCALE, HEIGHT/DOWNSCALE);

	let mut window = Window::new(
		"Test - ESC to exit",
		WIDTH,
		HEIGHT,
		WindowOptions::default(),
	)
		.unwrap_or_else(|e| {
			panic!("{}", e);
		});

	// Limit to max ~60 fps update rate
	window.limit_update_rate(Some(std::time::Duration::from_micros(16600)));

	while window.is_open() && !window.is_key_down(Key::Escape) {
		// Draw:
		let fluid:&FluidVolume = &fluid_tick;

		for (idx, pixel) in buffer.iter_mut().enumerate() {
			let x_screen = idx%WIDTH;
			let y_screen = idx/WIDTH;
			*pixel = (fluid.get_density(x_screen/DOWNSCALE, y_screen/DOWNSCALE) * (255f32)).min(255f32) as u32;
		}

		// Update Inputs and Simulate:
		let (mx, my) = window.get_mouse_pos(MouseMode::Clamp).unwrap();
		let btn = window.get_mouse_down(MouseButton::Left);
		if btn { fluid_tick.set_density(mx as usize/DOWNSCALE, my as usize/DOWNSCALE, INFLOW_RATE + fluid_tick.get_density(mx as usize, my as usize)); }

		for _ in 0..SUBSTEPS {
			fluid_tick.step_density_diffusion(DIFFUSION_RATE, SIMULATION_ITERATIONS, &mut fluid_tock);
			fluid_tock.step_density_with_velocity(1.0f32, &mut fluid_tick);
		}

		// We unwrap here as we want this code to exit if it fails. Real applications may want to handle this in a different way
		window
			.update_with_buffer(&buffer, WIDTH, HEIGHT)
			.unwrap();
	}
}
