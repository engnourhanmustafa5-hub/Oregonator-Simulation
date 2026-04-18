function lap = neumannLaplacian(U, dx)
    UPad = zeros(size(U,1)+2, size(U,2)+2);
    UPad(2:end-1, 2:end-1) = U;
    UPad(1,   2:end-1) = U(1,   :);
    UPad(end, 2:end-1) = U(end, :);
    UPad(2:end-1, 1)   = U(:,   1);
    UPad(2:end-1, end) = U(:,   end);
    UPad(1,   1)   = U(1,   1);
    UPad(1,   end) = U(1,   end);
    UPad(end, 1)   = U(end, 1);
    UPad(end, end) = U(end, end);

    lap = (UPad(2:end-1,3:end) + UPad(2:end-1,1:end-2) + ...
           UPad(3:end,2:end-1) + UPad(1:end-2,2:end-1) - ...
           4*UPad(2:end-1,2:end-1)) / (dx^2);
end
