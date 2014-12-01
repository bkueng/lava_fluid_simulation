

surface lava ( float Ka = 1, Kd = .5, Ks = .5, roughness = .1;
                         color specularcolor = 1;
						 string normal_tex = "";
						 string mask_tex = "";
						 string temperature_tex = "";
						 float min_temperature = 0.6;) {

	Ci = Cs;

	normal Nf = faceforward(normalize(N), I );
	point PNDC = transform("NDC", P);


	Oi = Os;
	//culling
	if (mask_tex != "") {
		float mask = float texture(mask_tex, xcomp(PNDC), ycomp(PNDC));
		float min_thres = 0.30;
		float max_thres = 0.99;
		if (mask < min_thres) {
			Ci = 0;
			Oi = 0; //set transparent
		} else {

			//non-zero: render surface
			if (normal_tex != "") {
				vector d = 0.5;
				Nf = normalize(vector texture(normal_tex, xcomp(PNDC), ycomp(PNDC))-d);

				//perturb randomly: problem: no temporal consistency
//				vector tangent, bitangent;
//				vector c1 = Nf ^ vector(0, 0, 1);
//				vector c2 = Nf ^ vector(0, 1, 0);
//				if (c1.c1 > c2.c2) {
//					tangent = c1;
//				} else {
//					tangent = c2;
//				}
//				tangent = normalize(tangent);
//				bitangent = normalize(Nf ^ tangent);
//				tangent += Nf * (float noise(P*400)-0.5)*0.3;
//				bitangent += Nf * (float noise(P*400)-0.5)*0.3;
//				Nf = normalize(tangent ^ bitangent);
//				Ks *= 1.2; /////////////////////////////////////////////////
			}
			float temperature = 1;
			if (temperature_tex != "") {
				//Color based on temperature
				temperature = float texture(temperature_tex, xcomp(PNDC), ycomp(PNDC));


				temperature = clamp((temperature-min_temperature)/(1-min_temperature), 0, 1);


				//add high freq noise
				//FIXME: this would look great but because the noise is based on
				//the position P, there is no temporal consistency. Thus for a
				//video is does not work! We would need some values that move
				//with the surface (temperature does not work)
//				temperature += (1-pow(temperature,5))*(noise(P*700)-0.5)*0.3;
//				temperature = min(1, abs(temperature));

				Ci = color spline(temperature,
						color(0.1294, 0.1137, 0.1059),/* lowest temperature */
						color(0.2000, 0.2000, 0.1725), /* dark gray */
						color(0.6667, 0.0706, 0.0118), /* more darker red */
						color(0.7686, 0.2000, 0.0078), /* darker red */
						color(0.9725, 0.1098, 0.0392), /* red */
						color(0.9725, 0.1098, 0.0392), /* red */
						color(0.9569, 0.5843, 0.0000), /* bright orange */
						color(0.9804, 0.9490, 0.2745), /* bright yellow */
						color(1, 0.996, 0.4),
						color(1, 0.996, 0.5),
						color(1, 1, 1) /* highest temperature */
						);
			}
			vector In = normalize(I);

			//make higher temperatures brighter
			color amb = ambient();
			float lerp_thres = 0.3;
			if (temperature > lerp_thres)
				amb = mix(amb, color(2), (temperature-lerp_thres)/(1-lerp_thres));

			Ci = Ci * (Ka*amb + Kd*diffuse(Nf))
				+ specularcolor * max(0, (1-temperature*3))*Ks*specular(Nf,-In,roughness);


		
			if (mask < max_thres) {
				Oi = (mask-min_thres)/(max_thres-min_thres);
				Ci *= Oi;
			}
		}
	}
}

