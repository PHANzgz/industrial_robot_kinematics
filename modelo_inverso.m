%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modelo geom�trico inverso para un robot de 6 ejes con mu�eca Roll-Pitch-Roll
%
% Q = MODELO_INVERSO(ROBOT, T) devuelve las coordenadas articulares que permiten
% a "robot" alcanzar la localizaci�n representada por la transformaci�n "T" (ya
% sea matriz o SE3). Representa una soluci�n anal�tica (enfoque algebr�ico o
% geom�trico), v�lida �nicamente para "robot".
%
% Q = MODELO_INVERSO(ROBOT, T, CONFIG_1,...) permite especificar en un n�
% variable de argumentos la configuraci�n codificada como una cadena de
% caracteres. Por ejemplo, para un robot angular t�pico (8 posibles configs):
%
% 'frente' (por defecto) o 'espaldas' para CONFIG_1
% 'arriba' (por defecto) o 'abajo' para CONFIG_2
% 'noflip' (por defecto) o 'flip' para CONFIG_3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function q = modelo_inverso(robot, T, varargin)

    %Definicion de un epsilon
    eps = 1e-6;
        
    % Si T no es de la clase SE3 (p ej una matriz 4x4), se convierte a SE3
    T = SE3(T);
 
    % Obtenci�n de los par�metros DH en la notaci�n usada en clase
    d = [robot.links.d];
    a = [robot.links.a];
    %ap = [robot.links.alpha];
    
    % Obtener 0T6, eliminando de la cadena cinem�tica las transformaciones
    % base y tool (estas transformaciones son por defecto la identidad)
    T = inv(robot.base) * T * inv(robot.tool);
    
    % Obtenci�n de los elementos de 0T6 en la notaci�n usada en clase
    tn = T.n; to = T.o; ta = T.a; p = T.t;
    
    % Tratamiento de las configuraciones
    %Config por defecto
    config_1 = 'frente';
    config_2 = 'arriba';
    config_3 = 'noflip';
    if nargin >= 3
        config_1 = lower(varargin{1});
    end
    if nargin >= 4
        config_2 = lower(varargin{2});
    end
    if nargin >= 5
        config_3 = lower(varargin{3});
    end
    
    % Determinaci�n de la posici�n del punto mu�eca
    m = p - d(6)*ta;
    
    q = zeros(1,6);
    % Comprobacion alcanzabilidad
    alc_xy = a(2) + sqrt(d(4)^2 + a(3)^2);
    alc_z = d(1) + a(2) + a(3) + d(4);
    
    if ( sqrt(m(1)^2 + m(2)^2) > alc_xy || m(3) > alc_z ) % salir si fuera de alcance
        disp("Posicion no alcanzable, devolviendo vector articular nulo");
        q(:) = 0;
        return;
    end
    
    % Soluci�n para q1
    if (contains(config_1,'espaldas'))
        disp("Configuracion de espaldas no disponible, devolviendo configuraci�n de frente.");
    end
    q(1) = atan2( m(2),m(1) );
    
    % Soluci�n para q2
    if (contains(config_2,'abajo'))
        disp("Configuracion codo abajo no disponible, devolviendo configuraci�n codo arriba.");
    end

    q2_atan2 = atan2(d(1) - m(3), (m(1))/cos(q(1)) - a(1));
    
    q2_acosnum = a(1)^2+a(2)^2-a(3)^2+d(1)^2-d(4)^2+m(3)^2-2*d(1)*m(3)+m(1)^2/ ...
        cos(q(1))^2-2*a(1)*m(1)/cos(q(1));
    
    q2_acosden = 2*a(2)*sqrt( (m(1)/cos(q(1)) - a(1))^2 +(d(1)- m(3))^2  );

    q(2) =  q2_atan2 - acos(q2_acosnum / q2_acosden );
    
    % Soluci�n para q3
    q3A = a(3)*cos(q(2)) - d(4)*sin(q(2));
    q3B = a(3)*sin(q(2)) + d(4)*cos(q(2));
    q3F = m(1)/cos(q(1)) - a(2)*cos(q(2)) - a(1);
    q3G = d(1) - a(2)*sin(q(2)) -m(3);
    q(3) = atan2( (q3A*q3G-q3B*q3F)/(q3A^2 + q3B^2), (q3A*q3F+q3B*q3G)/(q3A^2 + q3B^2));
    
    % Soluci�n para q5
    if (contains(config_3,'noflip'))
        q(5) = acos(-ta(3)*cos(q(2)+q(3)) - ta(1)*sin(q(2)+q(3))*cos(q(1)) - ...
            ta(2)*sin(q(2)+q(3))*sin(q(1)));
    else
        q(5) = -acos(-ta(3)*cos(q(2)+q(3)) - ta(1)*sin(q(2)+q(3))*cos(q(1)) - ...
            ta(2)*sin(q(2)+q(3))*sin(q(1)));
    end
    
    % Tratamiento singularidad mu�eca
    if ( abs(q(5)) < eps ) % alineaci�n
        disp(" Cerca o en singularidad de mu�eca, devolviendo una posible soluci�n.")
        q4q6 = atan2( tn(1)*sin(q(1)) - tn(2)*cos(q(1)), tn(1)*cos(q(2)+q(3))*...
                cos(q(1))-tn(3)*sin(q(2)+q(3))+tn(2)*cos(q(2)+q(3))*sin(q(1)) );
        % de las posibles soluciones, se da una arbitraria
        q(4) = 0;
        q(6) = q4q6;
        
    else % c�lculo normal
         % Soluci�n para q4
         q4num1 = -ta(1)*sin(q(1)) + ta(2)*cos(q(1));
         q4num2 = -ta(1)*cos(q(2)+q(3))*cos(q(1)) + ta(3)*sin(q(2)+q(3)) - ...
                ta(2)*cos(q(2)+q(3))*sin(q(1));
         q(4) = atan2( (q4num1/sin(q(5))), (q4num2/sin(q(5)) ) );
    
        % Soluci�n para q6
        q6num1 = to(3)*cos(q(2)+q(3))+to(1)*sin(q(2)+q(3))*cos(q(1)) + ...
                to(2)*sin(q(2)+q(3))*sin(q(1));
        q6num2 = -tn(3)*cos(q(2)+q(3))-tn(1)*sin(q(2)+q(3))*cos(q(1)) - ...
                tn(2)*sin(q(2)+q(3))*sin(q(1));
        q(6) = atan2( (q6num1/sin(q(5))), (q6num2/sin(q(5))));
    end
    
    % Tratamiento de articulaciones fuera de rango
    for link=1:6
        lim_i = robot.links(link).qlim;
        if ( q(link) < lim_i(1) || q(link) > lim_i(2) )
            disp("Articulaci�n " + string(link) + " fuera de rango; devolviendo valor 0 para todo el vector")
            q(:) = 0;
        end
    end
    
    % Formateo de la salida (redondear numeros a 0 si q<eps)
    q(abs(q)<eps) = 0;
    
end