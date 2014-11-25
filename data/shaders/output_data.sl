
surface output_data ( float do_mask=1; uniform float Temp = 0;
						output varying vector Nn = 0;
						 output varying float Mask = 0;
						 output varying float TempOut = 0;
						 output varying float Zdist = 0;) {
	//TODO: add temperature as output?
	// + particle identifier?

	TempOut = Temp;

	normal Nf = faceforward( normalize(N), I );
	Nn = normalize(Nf);
	point Pcam = transform("camera", P);
	Zdist = zcomp(Pcam);

	if (do_mask == 1)
		Mask = 1;
	else
		Mask = 0;

    Ci = Cs;
	Oi = 1;
}

