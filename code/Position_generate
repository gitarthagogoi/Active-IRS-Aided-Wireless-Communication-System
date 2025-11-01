function [Dis_BStoRIS, Dis_BStoUser, Dis_RIStoUser]=Position_generate(K,user_X)
Dis_BStoRIS=0;
Dis_BStoUser=zeros(K,1);
Dis_RIStoUser=zeros(K,1);

BS_position=[0 -60 0];
RIS_position=[300 10 0];

user_position=zeros(K,3);
% user_position(4,:)=[dist+cos(angle4) sin(angle4)];
user_position(1,:)=[user_X+5*cos(2*pi*rand()) 5*sin(2*pi*rand()) 0];
user_position(2,:)=[user_X+5*cos(2*pi*rand()) 5*sin(2*pi*rand()) 0];
user_position(3,:)=[user_X+5*cos(2*pi*rand()) 5*sin(2*pi*rand()) 0];
user_position(4,:)=[user_X+5*cos(2*pi*rand()) 5*sin(2*pi*rand()) 0];
% user_position(5,:)=[dist+4*(rand()-0.5) 4*(rand()-0.5)];
% user_position(6,:)=[dist+4*(rand()-0.5) 4*(rand()-0.5)];

Dis_BStoRIS=distance(BS_position,RIS_position);

for k=1:K
	user_position_temp=reshape(user_position(k,:),3,1);
	Dis_BStoUser(k)=distance(BS_position,user_position_temp);
end

for k=1:K
	user_position_temp=reshape(user_position(k,:),3,1);
	Dis_RIStoUser(k)=distance(RIS_position,user_position_temp);
end

end
