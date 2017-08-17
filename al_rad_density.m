function r_den = al_rad_density(zxy,r_edge)
    % YZ plane - radial density profile (X-integrated)
    r_cent=r_edge(1:end-1)+0.5*diff(r_edge);
    r_diff=diff(r_edge);
    
    r_norm_area=2*pi*r_cent.*r_diff;        % normalisation area

    % preprocess input
    yz=zxy(:,[3,1]);                 % cull X --> YZ
    
    % get radii
    r_vect=sqrt(sum(yz.^2,2));
    
    r_den=histcounts(r_vect,r_edge);     % get histogram counts
    r_den=r_den./r_norm_area;       % normalise counts by area    
end