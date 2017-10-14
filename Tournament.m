function w=Tournament(pop)

global nTournament;
h=randsample(numel(pop),nTournament);
[costorder,o]=min([pop(h).Cost]);
w=h(o);
end
 