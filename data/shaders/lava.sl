

surface lava ( float Ka = 1, Kd = .5, Ks = .5, roughness = .1;
                         color specularcolor = 1;
						 string normal_tex = "";
						 string mask_tex = "";
						 string temperature_tex = "";) {

	Ci = Cs;

	normal Nf = faceforward(normalize(N), I );
	point PNDC = transform("NDC", P);


	Oi = Os;
	//culling
	if (mask_tex != "") {
		float mask = float texture(mask_tex, xcomp(PNDC), ycomp(PNDC));
		float min_thres = 0.90;
		float max_thres = 0.99;
		if (mask < min_thres) {
			Ci = 0;
			Oi = 0; //set transparent
		} else {

			//non-zero: render surface
			if (normal_tex != "") {
				vector d = 0.5;
				Nf = normalize(vector texture(normal_tex, xcomp(PNDC), ycomp(PNDC))-d);
			}
			float temperature = 1;
			if (temperature_tex != "") {
				//Color based on temperature
				temperature = float texture(temperature_tex, xcomp(PNDC), ycomp(PNDC));

				//TODO: add high freq noise? -> must be temporal consistent!

				temperature = clamp(temperature*4-2.9, 0, 1);
				float temp = temperature;
				Ci = color spline(temp,
						color(0.1294, 0.1137, 0.1059),/* lowest temperature */
						color(0.2000, 0.2000, 0.1725), /* dark gray */
						color(0.6667, 0.0706, 0.0118), /* more darker red */
						color(0.7686, 0.2000, 0.0078), /* darker red */
						color(0.9725, 0.1098, 0.0392), /* red */
						color(0.9725, 0.1098, 0.0392), /* red */
						color(0.9569, 0.5843, 0.0000), /* bright orange */
						color(0.9804, 0.9490, 0.2745), /* bright yellow */
						color(1, 0.996, 0.8),
						color(1, 1, 1) /* highest temperature */
						);
			}
			vector In = normalize(I);

			//make higher temperatures brighter
			color amb = ambient();
			float lerp_thres = 0.3;
			if (temperature > lerp_thres)
				amb = mix(amb, color(1), (temperature-lerp_thres)/lerp_thres);

			Ci = Ci * (Ka*amb + Kd*diffuse(Nf))
				+ specularcolor * max(0, (1-temperature*3))*Ks*specular(Nf,-In,roughness);


		
			if (mask < max_thres) {
				Oi = (mask-min_thres)/(max_thres-min_thres);
				Ci *= Oi;
			}
		}
	}
}

