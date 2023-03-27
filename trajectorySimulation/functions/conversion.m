function output = conversion( input, from, to )
            au = 1.495978707E11;
            sy = seconds(years(1));
            g0 = 9.80665;
            
            switch from
                case "km"
                    switch to
                        case "au"
                            output = ( 1000 * input ) / au;
                        otherwise
                            error("unkown 'to' argument for 'from' argument.")
                    end

                case {'Isp','isp', 'specific impulse'}
                    switch to
                        case {'ve', 'Ve', 'exhaust velocity', 'effective exhaust velocity'}
                            output = input * g0;
                        otherwise
                            error("unkown 'to' argument for 'from' argument.")
                    end
                
                case {'ve', 'Ve', 'exhaust velocity', 'effective exhaust velocity'}
                    switch to
                        case {'Isp','isp', 'specific impulse'}
                            output = input / g0;
                        otherwise
                            error("unkown 'to' argument for 'from' argument.")
                    end

                case "m"
                    switch to
                        case "au"
                            output = ( input ) / au;
                        otherwise
                            error("unkown 'to' argument for 'from' argument.")
                    end
                    
               case "m^2"
                    switch to
                        case "au^2"
                            output = ( input ) / au^2;
                        otherwise
                            error("unkown 'to' argument for 'from' argument.")
                    end

               case "m^3"
                    switch to
                        case "au^3"
                            output = ( input ) / au^3;
                        otherwise
                            error("unkown 'to' argument for 'from' argument.")
                    end
                
               case "1/m^3"
                    switch to
                        case "1/au^3"
                            output = ( input ) * au^3;
                        otherwise
                            error("unkown 'to' argument for 'from' argument.")
                    end

                case "s"
                    switch to
                        case {"year","y","Y","Year"}
                            output = input / sy;
                        otherwise
                            error("unkown 'to' argument for 'from' argument.")
                    end

                case {"year","y","Y","Year"}
                    switch to
                        case {"s","second","seconds"}
                            output = input * sy;
                        otherwise
                            error("unkown 'to' argument for 'from' argument.")
                    end

                case "km/s"
                    switch to
                        case {"au/year","AU/Year","au/y"}
                            output = (1000 * input) / ( au / sy );
                        otherwise
                            error("unkown 'to' argument for 'from' argument.")
                    end

                case {"km/s^2","kilometer/second^2","kilometers/seconds^2"}
                    switch to
                        case {"au/year^2","AU/Year^2","au/y^2"}
                            output = (1000 * input) / ( au / sy^2 );
                        otherwise
                            error("unkown 'to' argument for 'from' argument.")
                    end

                case {"km^3/s^2","kilometer^3/second^2","kilometers^3/seconds^2"}
                    switch to
                        case {"au^3/year^2","AU^3/Year^2","au^3/y^2"}
                            output = (1000 * input) / ( au^3 / sy^2 );
                        otherwise
                            error("unkown 'to' argument for 'from' argument.")
                    end

                case "m/s"
                    switch to
                        case {"au/year","AU/Year","au/y"}
                            output = input / ( au / sy );
                        otherwise
                            error("unkown 'to' argument for 'from' argument.")
                    end

                case {"m/s^2","meter/second^2","meters/seconds^2"}
                    switch to
                        case {"au/year^2","AU/Year^2","au/y^2"}
                            output = input / ( au / sy^2 );
                        otherwise
                            error("unkown 'to' argument for 'from' argument.")
                    end

                case {"m^3/s^2","meter^3/second^2","meters^3/seconds^2"}
                    switch to
                        case {"au^3/year^2","AU^3/Year^2","au^3/y^2"}
                            output = input / ( au^3 / sy^2 );
                        otherwise
                            error("unkown 'to' argument for 'from' argument.")
                    end

                case {"au","AU","Astronomical Unit","Au"}
                    switch to
                        case {"m","meter","meters"}
                            output = input * au;
                        case {"km","kilometer","kilometers"}
                            output = ( input * au ) / 1000;
                        otherwise
                            error("unkown 'to' argument for 'from' argument.")
                    end

                case {"au/y","AU/Year","AU/y","au/Year","AU/year","au/year"}
                    switch to
                        case {"m/s","meter/second","meters/seconds"}
                            output = input * ( au / sy );
                        case {"km/s","kilometer/second","kilometers/seconds"}
                            output = ( input * ( au / sy ) ) / 1000;
                        otherwise
                            error("unkown 'to' argument for 'from' argument.")
                    end

                case {"au/y^2","AU/Year^2","AU/y^2","au/Year^2","AU/year^2","au/year^2"}
                    switch to
                        case {"m/s^2","meter/second^2","meters/seconds^2"}
                            output = input * ( au / sy^2 );
                        case {"km/s^2","kilometer/second^2","kilometers/seconds^2"}
                            output = ( input * ( au / sy^2 ) ) / 1000;
                        otherwise
                            error("unkown 'to' argument for 'from' argument.")
                    end

                otherwise
                    error(" unknown 'from' argument.")
            end
        end