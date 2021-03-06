JEOS = @(rho,T) rho.*T+rho.^2.*(JohnsonParams.x1.*T+JohnsonParams.x2.*sqrt(T)+JohnsonParams.x3+JohnsonParams.x4./T+JohnsonParams.x5./T.^2) ...
             + rho.^3.*(JohnsonParams.x6.*T+JohnsonParams.x7+JohnsonParams.x8./T+JohnsonParams.x9./T.^2) + rho.^4.*(JohnsonParams.x10.*T+JohnsonParams.x11+JohnsonParams.x12./T) ...
             + rho.^5.*(JohnsonParams.x13) + rho.^6.*(JohnsonParams.x14./T+JohnsonParams.x15./T.^2) + rho.^7.*(JohnsonParams.x16./T) ...
             + rho.^8.*(JohnsonParams.x17./T+JohnsonParams.x18./T.^2) + rho.^9.*(JohnsonParams.x19./T.^2) ...
             + ( rho.^3.*(JohnsonParams.x20./T.^2+JohnsonParams.x21./T.^3) ...
               + rho.^5.*(JohnsonParams.x22./T.^2+JohnsonParams.x23./T.^4)   + rho.^7.*(JohnsonParams.x24./T.^2+JohnsonParams.x25./T.^3) ...
               + rho.^9.*(JohnsonParams.x26./T.^2+JohnsonParams.x27./T.^4)   + rho.^11.*(JohnsonParams.x28./T.^2+JohnsonParams.x29./T.^3) ...
               + rho.^13.*(JohnsonParams.x30./T.^2+JohnsonParams.x31./T.^3+JohnsonParams.x32./T.^4) ) .* exp(-JohnsonParams.g.*rho.^2);