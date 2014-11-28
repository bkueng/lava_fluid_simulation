
surface paintedplasticglow ( float Ka = 1, Kd = .5, Ks = .5, roughness = .1;
                         color specularcolor = 1;
                         string texturename = "";
						 string glow_tex = "";
						 float glow_amp = 1;) {
    Ci = Cs;
    if (texturename != "")
		Ci *= color texture (texturename);

    normal Nf = faceforward (normalize(N),I);
    Ci = Ci * (Ka*ambient() + Kd*diffuse(Nf)) + specularcolor * Ks*specular(Nf,-normalize(I),roughness);
	
	if (glow_tex != "") {
		point PNDC = transform("NDC", P);
		color glow_color = color texture(glow_tex, xcomp(PNDC), ycomp(PNDC));
		Ci += clamp(glow_amp * glow_color, 0, 1);
	}
    Oi = Os;
	Ci *= Oi;
}

