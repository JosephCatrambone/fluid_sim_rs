
mod fluid;

use minifb::{Key, MouseButton, MouseMode, Window, WindowOptions};
use crate::fluid::FluidVolume;

const WIDTH: usize = 640;
const HEIGHT: usize = 640;
const DOWNSCALE: usize = 16;
const DIFFUSION_RATE: f32 = 1.0f32;
const SIMULATION_ITERATIONS: usize = 10;
const INFLOW_RATE: f32 = 10.0f32;

fn main() {
	let mut buffer: Vec<u32> = vec![0; WIDTH * HEIGHT];
	let mut fluid = fluid::FluidVolume::new(WIDTH/DOWNSCALE, HEIGHT/DOWNSCALE);
	//let mut last_mouse: (f32, f32) = (0.0, 0.0);

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
		let left_btn = window.get_mouse_down(MouseButton::Left);
		if left_btn { fluid.set_density(mx as usize/DOWNSCALE, my as usize/DOWNSCALE, INFLOW_RATE + fluid.get_density(mx as usize, my as usize)); }
		let right_btn = window.get_mouse_down(MouseButton::Right);
		if right_btn {
			//let mdx = mx - last_mouse.0;
			//let mdy = my - last_mouse.1;
			//last_mouse = (mx, my);
			let mdx = mx - (WIDTH/2) as f32;
			let mdy = my - (HEIGHT/2) as f32;

			println!("Fluid velocity: {} {}", &mdx, &mdy);
			fluid.set_velocity((WIDTH/2)/DOWNSCALE, (HEIGHT/2)/DOWNSCALE, (mdx, mdy));
		}

		fluid.step(DIFFUSION_RATE, 0.1f32, SIMULATION_ITERATIONS);

		// We unwrap here as we want this code to exit if it fails. Real applications may want to handle this in a different way
		window
			.update_with_buffer(&buffer, WIDTH, HEIGHT)
			.unwrap();
	}
}
