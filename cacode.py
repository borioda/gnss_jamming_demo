def cacode( svid ) :
    import numpy as n
    # Summary:
    #   function that generates the code associated to a specific PRN
    #   in this case the GPS L1 C/A code is considered
    #
    # Arguments:
    #   svid:   satellite identifier
    #
    # Returns:
    #   code - vector of 1023 elements

    # Taps assigment for the different svn
    x = n.array([[2, 6],[3, 7],[4, 8],[5, 9],[1, 9],\
                [2, 10],[1, 8],[2, 9],[3, 10],[2, 3],\
                [3, 4],[5, 6],[6, 7],[7, 8],[8, 9],\
                [9, 10],[1, 4],[2, 5],[3, 6],[4, 7],\
                [5, 8],[6, 9],[1, 3],[4, 6],[5, 7],\
                [6, 8],[7, 9],[8, 10],[1, 6],[2, 7],\
                [3, 8],[4, 9]])
            
    # select the taps corresponding to the selected SV
    tap1 = x[ svid-1, 0]-1       # remember python indexing starts from 0!
    tap2 = x[ svid-1, 1]-1

    # generate the C/A cpde
    g1 = n.ones(10)
    g2 = n.ones(10)

    code = n.zeros(1023)
    
    for ii in range(0, 1023, 1) :
        code[ ii ] = n.logical_xor( n.logical_xor( g1[ 9 ], g2[ tap1 ] ), g2[ tap2] )
        temp = n.logical_xor( g1[ 2 ], g1[ 9 ] )      
        g1[ 1:10 ] = g1[ 0:9 ]
        g1[ 0 ] = n.double( temp )
        temp = n.logical_xor( g2[1],n.logical_xor(g2[2],n.logical_xor(g2[5],\
               n.logical_xor(g2[7],n.logical_xor(g2[8],g2[9])))))
        g2[ 1:10 ] = g2[ 0:9 ]
        g2[ 0 ] = n.double( temp )
    # end for
    
    # make the sequence bipolar
    code[ : ] = 1 - 2 * code[:]
    return code         
           

