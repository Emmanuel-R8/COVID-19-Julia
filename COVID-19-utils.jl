#-------------------------------------------------------------------------------------------------
#--
#-- Linear interpolation of Y values for a new set of X values.
#--
function linearInterpolation(x, xvar, yvar)
    l = length(xvar)

    # If l = 1, ratio will be the only one
    if l == 1
        return yvar[1]
    else
        for i in 2:l
            x1 = xvar[i-1]
            x2 = xvar[i  ]

            if x < x2
                deltaY = yvar[i] - yvar[i-1]
                return yvar[i-1] + deltaY * (x - x1) / (x2 - x1)
            end
        end

        # Last possible choice
        return yvar[l]
    end
end
