

surface lava ( float Ka = 1, Kd = .5, Ks = .5, roughness = .1;
                         color specularcolor = 1;
						 string normal_tex = "";
						 string mask_tex = "";
						 string temperature_tex = "";) {
	//TODO: particle identifier texture?

	Ci = Cs;

	normal Nf = faceforward(normalize(N), I );
	point PNDC = transform("NDC", P);

	if (normal_tex != "") {
		vector d = 0.5;
		Nf = normalize(vector texture(normal_tex, xcomp(PNDC), ycomp(PNDC))-d);
	}
	vector In = normalize(I);
	Ci = Ci * (Ka*ambient() + Kd*diffuse(Nf)) + specularcolor * Ks*specular(Nf,-In,roughness);


	Oi = Os;
	//culling
	if (mask_tex != "") {
		float mask = float texture(mask_tex, xcomp(PNDC), ycomp(PNDC));
		if (mask < 0.95) {
			Ci = 0;
			Oi = 0; //set transparent
		}
	}
}

