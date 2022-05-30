
mod fluid;

use minifb::{Key, MouseButton, MouseMode, Window, WindowOptions};
use crate::fluid::FluidVolume;

const WIDTH: usize = 640;
const HEIGHT: usize = 640;
const DOWNSCALE: usize = 16;
const DIFFUSION_RATE: f32 = 0.5f32;
const SIMULATION_ITERATIONS: usize = 10;
const INFLOW_RATE: f32 = 10.0f32;

fn main() {
	let mut buffer: Vec<u32> = vec![0; WIDTH * HEIGHT];
	let mut fluid = fluid::FluidVolume::new(WIDTH/DOWNSCALE, HEIGHT/DOWNSCALE);
	let mut last_mouse: (f32, f32) = (0.0, 0.0);

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
		for (idx, pixel) in buffer.iter_mut().enumerate() {
			let x_screen = idx%WIDTH;
			let y_screen = idx/WIDTH;
			let pixel_blue = ((fluid.get_density(x_screen/DOWNSCALE, y_screen/DOWNSCALE) * (255f32))).min(255f32) as u32;
			let pixel_red = ((fluid.get_velocity_x(x_screen/DOWNSCALE, y_screen/DOWNSCALE) * (255f32))).abs().min(255f32) as u32;
			let pixel_green = ((fluid.get_velocity_y(x_screen/DOWNSCALE, y_screen/DOWNSCALE) * (255f32))).abs().min(255f32) as u32;
			*pixel = ((0x010000 * pixel_red) + (0x000100 * pixel_green) + (0x000001 * pixel_blue));
		}

		// Update Inputs and Simulate:
		let (mx, my) = window.get_mouse_pos(MouseMode::Clamp).unwrap();
		let mdx = mx - last_mouse.0;
		let mdy = my - last_mouse.1;
		last_mouse = (mx, my);

		if window.get_mouse_down(MouseButton::Left) {
			fluid.set_density(mx as usize/DOWNSCALE, my as usize/DOWNSCALE, INFLOW_RATE + fluid.get_density(mx as usize/DOWNSCALE, my as usize/DOWNSCALE));
		}
		if window.get_mouse_down(MouseButton::Right) {
			fluid.set_velocity(mx as usize/DOWNSCALE, my as usize/DOWNSCALE, (mdx*0.01f32, mdy*0.01f32));
		}

		fluid.step(DIFFUSION_RATE, 0.1f32, SIMULATION_ITERATIONS);

		// We unwrap here as we want this code to exit if it fails. Real applications may want to handle this in a different way
		window
			.update_with_buffer(&buffer, WIDTH, HEIGHT)
			.unwrap();
	}
}
