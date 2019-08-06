function cdata = face2vertex(cdata,faces,nvert)

    fmax = max(faces(:));
    if nargin < 3, nvert=fmax; end
    if size(faces,1)~=3, faces=faces'; end

    assert( size(faces,1)==3, 'Bad faces size.' );
    assert( size(faces,2)==numel(cdata), 'Input size mismatch.' );
    assert( nvert >= fmax, 'Number of vertices too small.' );

    faces = faces(:);
    cdata = repelem( cdata(:), 3 ); % triplicate face colors

    nfpv  = accumarray( faces, 1, [nvert,1] ); % #of faces per vertex
    cdata = accumarray( faces, cdata, [nvert,1] ) ./ max(1,nfpv);

end