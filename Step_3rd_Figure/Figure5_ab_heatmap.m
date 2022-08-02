clear

load('\Weight_sig_Comp1.mat');
net_roi = importdata("\net_roi.csv")
network_name = net_roi.textdata(2:end,1);
network_name1 = unique(network_name);
network_name1(12) = [];
network_name1(14)={'Subcortical'};
network_name_new = {'Aud','CO','CP','DMN','DAN','FP','RSP','SM','SMla','Sal','VAN','Vis','Sub'};
[a,b] = find(Brain_roi_weight_sig<0.05);
Brain_roi_weight05 = zeros(352,352);
for i = 1:length(a)
   Brain_roi_weight05(a(i),b(i)) = fc_roi_MedianWeight(a(i),b(i));
end 
for i = 1:length(network_name1)
   index1 = find(strcmp(network_name,network_name1(i))==1);
   indexa(i) = index1(end);
   indexb(i)= index1(1);
   index(i) = round((indexa(i)+indexb(i))/2);
end
indexac =indexa;
indexac(9)=[];
indexbc = indexb;
indexbc(10)=[];
indexc = round((indexac+indexbc)/2);
[r,c] = find(Brain_net_weight_sig<0.05);
Brain_net_weight05 = zeros(14,14);
for i = 1:length(r)
    Brain_net_weight05_abs(r(i),c(i)) = fc_net_MedianWeight(r(i),c(i));
end
T = table(Brain_net_weight05_abs)
writetable(T, 'D:\projects\ABCD\EF\PLS\results\ef_task23\combat\2fold\comp1\net_abs.txt')
end
Brain_roi_weight05n = zeros(352,352);
for i = 1:length(r)
    startr = indexb(r(i));
    endr = indexa(r(i));
    row = [startr:endr];
    startc = indexb(c(i));
    endc = indexa(c(i));
    col = [startc:endc];
    Brain_roi_weight05n(row,col)= fc_roi_MedianWeight(row,col);
end

map = [0.537254929542542 0.909803926944733 0.901960790157318;0.535000026226044 0.906176447868347 0.903039216995239;0.532745122909546 0.902549028396606 0.90411764383316;0.530490219593048 0.898921549320221 0.905196070671082;0.52823531627655 0.89529412984848 0.906274497509003;0.525980412960052 0.891666650772095 0.907352924346924;0.523725509643555 0.888039231300354 0.908431351184845;0.521470606327057 0.884411752223969 0.909509837627411;0.519215703010559 0.880784332752228 0.910588264465332;0.516960799694061 0.877156853675842 0.911666691303253;0.514705896377563 0.873529434204102 0.912745118141174;0.512450993061066 0.869901955127716 0.913823544979095;0.510196089744568 0.866274535655975 0.914901971817017;0.50794118642807 0.86264705657959 0.915980398654938;0.505686283111572 0.859019637107849 0.917058825492859;0.503431379795074 0.855392158031464 0.91813725233078;0.501176476478577 0.851764738559723 0.919215679168701;0.498921602964401 0.848137259483337 0.920294106006622;0.496666699647903 0.844509840011597 0.921372532844543;0.494411796331406 0.840882360935211 0.922450959682465;0.492156893014908 0.837254881858826 0.923529386520386;0.48990198969841 0.833627462387085 0.924607872962952;0.487647086381912 0.829999983310699 0.925686299800873;0.485392183065414 0.826372563838959 0.926764726638794;0.483137279748917 0.822745084762573 0.927843153476715;0.480882376432419 0.819117665290833 0.928921580314636;0.478627473115921 0.815490186214447 0.930000007152557;0.476372569799423 0.811862766742706 0.931078433990479;0.474117666482925 0.808235287666321 0.9321568608284;0.471862763166428 0.80460786819458 0.933235287666321;0.46960785984993 0.800980389118195 0.934313714504242;0.467352956533432 0.797352969646454 0.935392141342163;0.465098053216934 0.793725490570068 0.936470568180084;0.462843149900436 0.790098071098328 0.937548995018005;0.460588246583939 0.786470592021942 0.938627481460571;0.458333343267441 0.782843172550201 0.939705908298492;0.456078439950943 0.779215693473816 0.940784335136414;0.453823536634445 0.775588274002075 0.941862761974335;0.451568633317947 0.77196079492569 0.942941188812256;0.44931373000145 0.768333375453949 0.944019615650177;0.447058856487274 0.764705896377563 0.945098042488098;0.444803953170776 0.761078417301178 0.946176469326019;0.442549049854279 0.757450997829437 0.94725489616394;0.440294146537781 0.753823518753052 0.948333323001862;0.438039243221283 0.750196099281311 0.949411749839783;0.435784339904785 0.746568620204926 0.950490176677704;0.433529436588287 0.742941200733185 0.951568603515625;0.43127453327179 0.739313721656799 0.952647089958191;0.429019629955292 0.735686302185059 0.953725516796112;0.426764726638794 0.732058823108673 0.954803943634033;0.424509823322296 0.728431403636932 0.955882370471954;0.422254920005798 0.724803924560547 0.956960797309875;0.420000016689301 0.721176505088806 0.958039224147797;0.417745113372803 0.717549026012421 0.959117650985718;0.415490210056305 0.71392160654068 0.960196077823639;0.413235306739807 0.710294127464294 0.96127450466156;0.410980403423309 0.706666707992554 0.962352931499481;0.408725500106812 0.703039228916168 0.963431358337402;0.406470596790314 0.699411809444427 0.964509785175323;0.404215693473816 0.695784330368042 0.965588212013245;0.401960790157318 0.692156910896301 0.966666698455811;0.39970588684082 0.688529431819916 0.967745125293732;0.397450983524323 0.68490195274353 0.968823552131653;0.395196080207825 0.68127453327179 0.969901978969574;0.392941176891327 0.677647054195404 0.970980405807495;0.390686273574829 0.674019634723663 0.972058832645416;0.388431370258331 0.670392155647278 0.973137259483337;0.386176496744156 0.666764736175537 0.974215686321259;0.383921593427658 0.663137257099152 0.97529411315918;0.38166669011116 0.659509837627411 0.976372539997101;0.379411786794662 0.655882358551025 0.977450966835022;0.377156883478165 0.652254939079285 0.978529393672943;0.374901980161667 0.648627460002899 0.979607820510864;0.372647076845169 0.645000040531158 0.980686247348785;0.370392173528671 0.641372561454773 0.981764733791351;0.368137270212173 0.637745141983032 0.982843160629272;0.365882366895676 0.634117662906647 0.983921587467194;0.363627463579178 0.630490243434906 0.985000014305115;0.36137256026268 0.626862764358521 0.986078441143036;0.359117656946182 0.62323534488678 0.987156867980957;0.356862753629684 0.619607865810394 0.988235294818878;0.317211329936981 0.550762534141541 0.989542484283447;0.277559906244278 0.481917232275009 0.990849673748016;0.237908497452736 0.413071900606155 0.992156863212585;0.198257088661194 0.344226598739624 0.993464052677155;0.158605664968491 0.27538126707077 0.994771242141724;0.118954248726368 0.206535950303078 0.996078431606293;0.0793028324842453 0.137690633535385 0.997385621070862;0.0396514162421227 0.0688453167676926 0.998692810535431;0 0 1;0 0 0.922291040420532;0 0 0.844582140445709;0 0 0.766873240470886;0 0 0.689164280891418;0 0 0.611455321311951;0 0 0.533746421337128;0 0 0.456037491559982;0 0 0.378328561782837;0 0 0.300619632005692;0 0 0.297387838363647;0 0 0.294156044721603;0 0 0.290924280881882;0 0 0.287692487239838;0 0 0.284460693597794;0 0 0.28122889995575;0 0 0.277997136116028;0 0 0.274765342473984;0 0 0.27153354883194;0 0 0.268301755189896;0 0 0.265069991350174;0 0 0.26183819770813;0 0 0.258606404066086;0 0 0.255374610424042;0 0 0.25214284658432;0 0 0.248911052942276;0 0 0.245679259300232;0 0 0.242447480559349;0 0 0.239215686917305;0.00219298247247934 0 0.226625382900238;0.00438596494495869 0 0.214035093784332;0.00657894741743803 0 0.201444789767265;0.00877192988991737 0 0.188854485750198;0.0109649123623967 0 0.176264196634293;0.0131578948348761 0 0.163673892617226;0.0153508773073554 0 0.151083588600159;0.0175438597798347 0 0.138493299484253;0.0197368431836367 0 0.125902995467186;0.0219298247247934 0 0.113312691450119;0.0241228081285954 0 0.100722394883633;0.0263157896697521 0 0.0881320983171463;0.028508773073554 0 0.0755417943000793;0.0307017546147108 0 0.062951497733593;0.0328947380185127 0 0.0503611974418163;0.0350877195596695 0 0.0377708971500397;0.0372807011008263 0 0.0251805987209082;0.0394736863672733 0 0.0125902993604541;0.0416666679084301 0 0;0.0949074104428291 0 0;0.148148149251938 0 0;0.201388895511627 0 0;0.254629641771317 0 0;0.307870358228683 0 0;0.361111104488373 0 0;0.414351850748062 0 0;0.467592597007751 0 0;0.520833313465118 0 0;0.57407408952713 0 0;0.627314805984497 0 0;0.680555582046509 0 0;0.733796298503876 0 0;0.787037014961243 0 0;0.840277791023254 0 0;0.893518507480621 0 0;0.946759283542633 0 0;1 0 0;0.999411761760712 0.0407843142747879 0.00941176526248455;0.998823523521423 0.0815686285495758 0.0188235305249691;0.998235285282135 0.122352942824364 0.0282352939248085;0.997647047042847 0.163137257099152 0.0376470610499382;0.997058808803558 0.20392157137394 0.0470588244497776;0.99647057056427 0.244705885648727 0.056470587849617;0.995882332324982 0.285490214824677 0.0658823549747467;0.995294094085693 0.326274514198303 0.0752941220998764;0.994705855846405 0.36705881357193 0.0847058817744255;0.994117617607117 0.407843142747879 0.0941176488995552;0.993529438972473 0.448627471923828 0.103529416024685;0.992941200733185 0.489411771297455 0.112941175699234;0.992352962493896 0.530196070671082 0.122352942824364;0.991764724254608 0.570980429649353 0.131764709949493;0.99117648601532 0.61176472902298 0.141176477074623;0.990588247776031 0.652549028396606 0.150588244199753;0.990000009536743 0.693333327770233 0.159999996423721;0.989411771297455 0.73411762714386 0.169411763548851;0.988823533058167 0.774901986122131 0.178823530673981;0.988235294818878 0.815686285495758 0.18823529779911;0.988380551338196 0.817961752414703 0.195303797721863;0.988525807857513 0.820237219333649 0.202372312545776;0.988671004772186 0.822512745857239 0.209440812468529;0.988816261291504 0.824788212776184 0.216509327292442;0.988961517810822 0.827063679695129 0.223577827215195;0.989106774330139 0.829339146614075 0.230646342039108;0.989251971244812 0.83161461353302 0.237714841961861;0.98939722776413 0.833890080451965 0.244783356785774;0.989542484283447 0.836165606975555 0.251851856708527;0.989687740802765 0.838441073894501 0.25892037153244;0.989832997322083 0.840716540813446 0.265988856554031;0.989978194236755 0.842992007732391 0.273057371377945;0.990123450756073 0.845267474651337 0.280125886201859;0.990268707275391 0.847543001174927 0.287194401025772;0.990413963794708 0.849818468093872 0.294262886047363;0.990559160709381 0.852093935012817 0.301331400871277;0.990704417228699 0.854369401931763 0.30839991569519;0.990849673748016 0.856644868850708 0.315468430519104;0.990994930267334 0.858920395374298 0.322536915540695;0.991140186786652 0.861195862293243 0.329605430364609;0.991285383701324 0.863471329212189 0.336673945188522;0.991430640220642 0.865746796131134 0.343742430210114;0.99157589673996 0.868022263050079 0.350810945034027;0.991721153259277 0.870297729969025 0.357879459857941;0.99186635017395 0.872573256492615 0.364947974681854;0.992011606693268 0.87484872341156 0.372016459703445;0.992156863212585 0.877124190330505 0.379084974527359;0.992302119731903 0.879399657249451 0.386153489351273;0.992447376251221 0.881675124168396 0.393221974372864;0.992592573165894 0.883950650691986 0.400290489196777;0.992737829685211 0.886226117610931 0.407359004020691;0.992883086204529 0.888501584529877 0.414427518844604;0.993028342723846 0.890777051448822 0.421496003866196;0.993173539638519 0.893052518367767 0.428564518690109;0.993318796157837 0.895327985286713 0.435633033514023;0.993464052677155 0.897603511810303 0.442701518535614;0.993609309196472 0.899878978729248 0.449770033359528;0.99375456571579 0.902154445648193 0.456838548183441;0.993899762630463 0.904429912567139 0.463907063007355;0.99404501914978 0.906705379486084 0.470975548028946;0.994190275669098 0.908980906009674 0.478044062852859;0.994335532188416 0.911256372928619 0.485112577676773;0.994480729103088 0.913531839847565 0.492181092500687;0.994625985622406 0.91580730676651 0.499249577522278;0.994771242141724 0.918082773685455 0.506318092346191;0.994916498661041 0.920358300209045 0.513386607170105;0.995061755180359 0.922633767127991 0.520455121994019;0.995206952095032 0.924909234046936 0.527523636817932;0.995352208614349 0.927184700965881 0.534592092037201;0.995497465133667 0.929460167884827 0.541660606861115;0.995642721652985 0.931735634803772 0.548729121685028;0.995787918567657 0.934011161327362 0.555797636508942;0.995933175086975 0.936286628246307 0.562866151332855;0.996078431606293 0.938562095165253 0.569934666156769;0.99622368812561 0.940837562084198 0.577003180980682;0.996368944644928 0.943113029003143 0.584071636199951;0.996514141559601 0.945388555526733 0.591140151023865;0.996659398078918 0.947664022445679 0.598208665847778;0.996804654598236 0.949939489364624 0.605277180671692;0.996949911117554 0.952214956283569 0.612345695495605;0.997095108032227 0.954490423202515 0.619414210319519;0.997240364551544 0.95676589012146 0.626482725143433;0.997385621070862 0.95904141664505 0.633551239967346;0.997530877590179 0.961316883563995 0.640619695186615;0.997676134109497 0.963592350482941 0.647688210010529;0.99782133102417 0.965867817401886 0.654756724834442;0.997966587543488 0.968143284320831 0.661825239658356;0.998111844062805 0.970418810844421 0.668893754482269;0.998257100582123 0.972694277763367 0.675962269306183;0.998402297496796 0.974969744682312 0.683030784130096;0.998547554016113 0.977245211601257 0.690099239349365;0.998692810535431 0.979520678520203 0.697167754173279;0.998838067054749 0.981796205043793 0.704236268997192;0.998983323574066 0.984071671962738 0.711304783821106;0.999128520488739 0.986347138881683 0.71837329864502;0.999273777008057 0.988622605800629 0.725441813468933;0.999419033527374 0.990898072719574 0.732510328292847;0.999564290046692 0.993173539638519 0.739578783512115;0.999709486961365 0.995449066162109 0.746647298336029;0.999854743480682 0.997724533081055 0.753715813159943;1 1 0.760784327983856];  
map05 = [0.537254929542542 0.909803926944733 0.901960790157318;0.535000026226044 0.906176447868347 0.903039216995239;0.532745122909546 0.902549028396606 0.90411764383316;0.530490219593048 0.898921549320221 0.905196070671082;0.52823531627655 0.89529412984848 0.906274497509003;0.525980412960052 0.891666650772095 0.907352924346924;0.523725509643555 0.888039231300354 0.908431351184845;0.521470606327057 0.884411752223969 0.909509837627411;0.519215703010559 0.880784332752228 0.910588264465332;0.516960799694061 0.877156853675842 0.911666691303253;0.514705896377563 0.873529434204102 0.912745118141174;0.512450993061066 0.869901955127716 0.913823544979095;0.510196089744568 0.866274535655975 0.914901971817017;0.50794118642807 0.86264705657959 0.915980398654938;0.505686283111572 0.859019637107849 0.917058825492859;0.503431379795074 0.855392158031464 0.91813725233078;0.501176476478577 0.851764738559723 0.919215679168701;0.498921602964401 0.848137259483337 0.920294106006622;0.496666699647903 0.844509840011597 0.921372532844543;0.494411796331406 0.840882360935211 0.922450959682465;0.492156893014908 0.837254881858826 0.923529386520386;0.48990198969841 0.833627462387085 0.924607872962952;0.487647086381912 0.829999983310699 0.925686299800873;0.485392183065414 0.826372563838959 0.926764726638794;0.483137279748917 0.822745084762573 0.927843153476715;0.480882376432419 0.819117665290833 0.928921580314636;0.478627473115921 0.815490186214447 0.930000007152557;0.476372569799423 0.811862766742706 0.931078433990479;0.474117666482925 0.808235287666321 0.9321568608284;0.471862763166428 0.80460786819458 0.933235287666321;0.46960785984993 0.800980389118195 0.934313714504242;0.467352956533432 0.797352969646454 0.935392141342163;0.465098053216934 0.793725490570068 0.936470568180084;0.462843149900436 0.790098071098328 0.937548995018005;0.460588246583939 0.786470592021942 0.938627481460571;0.458333343267441 0.782843172550201 0.939705908298492;0.456078439950943 0.779215693473816 0.940784335136414;0.453823536634445 0.775588274002075 0.941862761974335;0.451568633317947 0.77196079492569 0.942941188812256;0.44931373000145 0.768333375453949 0.944019615650177;0.447058856487274 0.764705896377563 0.945098042488098;0.444803953170776 0.761078417301178 0.946176469326019;0.442549049854279 0.757450997829437 0.94725489616394;0.440294146537781 0.753823518753052 0.948333323001862;0.438039243221283 0.750196099281311 0.949411749839783;0.435784339904785 0.746568620204926 0.950490176677704;0.433529436588287 0.742941200733185 0.951568603515625;0.43127453327179 0.739313721656799 0.952647089958191;0.429019629955292 0.735686302185059 0.953725516796112;0.426764726638794 0.732058823108673 0.954803943634033;0.424509823322296 0.728431403636932 0.955882370471954;0.422254920005798 0.724803924560547 0.956960797309875;0.420000016689301 0.721176505088806 0.958039224147797;0.417745113372803 0.717549026012421 0.959117650985718;0.415490210056305 0.71392160654068 0.960196077823639;0.413235306739807 0.710294127464294 0.96127450466156;0.410980403423309 0.706666707992554 0.962352931499481;0.408725500106812 0.703039228916168 0.963431358337402;0.406470596790314 0.699411809444427 0.964509785175323;0.404215693473816 0.695784330368042 0.965588212013245;0.401960790157318 0.692156910896301 0.966666698455811;0.39970588684082 0.688529431819916 0.967745125293732;0.397450983524323 0.68490195274353 0.968823552131653;0.395196080207825 0.68127453327179 0.969901978969574;0.392941176891327 0.677647054195404 0.970980405807495;0.390686273574829 0.674019634723663 0.972058832645416;0.388431370258331 0.670392155647278 0.973137259483337;0.386176496744156 0.666764736175537 0.974215686321259;0.383921593427658 0.663137257099152 0.97529411315918;0.38166669011116 0.659509837627411 0.976372539997101;0.379411786794662 0.655882358551025 0.977450966835022;0.377156883478165 0.652254939079285 0.978529393672943;0.374901980161667 0.648627460002899 0.979607820510864;0.372647076845169 0.645000040531158 0.980686247348785;0.370392173528671 0.641372561454773 0.981764733791351;0.368137270212173 0.637745141983032 0.982843160629272;0.365882366895676 0.634117662906647 0.983921587467194;0.363627463579178 0.630490243434906 0.985000014305115;0.36137256026268 0.626862764358521 0.986078441143036;0.359117656946182 0.62323534488678 0.987156867980957;0.356862753629684 0.619607865810394 0.988235294818878;0.337037056684494 0.58518522977829 0.98888885974884;0.317211329936981 0.550762534141541 0.989542484283447;0.297385632991791 0.516339898109436 0.990196108818054;0.277559906244278 0.481917232275009 0.990849673748016;0.257734209299088 0.447494566440582 0.991503238677979;0.237908497452736 0.413071900606155 0.992156863212585;0.218082800507545 0.378649264574051 0.992810487747192;0.198257088661194 0.344226598739624 0.993464052677155;0.178431376814842 0.309803932905197 0.994117617607117;0.158605664968491 0.27538126707077 0.994771242141724;0.138779953122139 0.240958616137505 0.995424866676331;0.118954248726368 0.206535950303078 0.996078431606293;0.0991285443305969 0.172113299369812 0.996731996536255;0.0793028324842453 0.137690633535385 0.997385621070862;0.059477124363184 0.103267975151539 0.998039245605469;0.0396514162421227 0.0688453167676926 0.998692810535431;0.0198257081210613 0.0344226583838463 0.999346375465393;0 0 1;0 0 0.990711331367493;0 0 0.98142272233963;0 0 0.972134113311768;0 0 0.96284544467926;0 0 0.861997306346893;0 0 0.761149108409882;0 0 0.660300970077515;0 0 0.559452831745148;0 0 0.458604663610458;0 0 0.357756525278091;0 0 0.256908357143402;0 0 0.256139099597931;0 0 0.255369871854782;0 0 0.254600614309311;0 0 0.25383135676384;0 0 0.253062129020691;0 0 0.25229287147522;0 0 0.251523643732071;0 0 0.2507543861866;0 0 0.24998514354229;0 0 0.249215885996819;0 0 0.248446643352509;0 0 0.247677400708199;0 0 0.246908158063889;0 0 0.246138900518417;0 0 0.245369657874107;0 0 0.244600415229797;0 0 0.243831172585487;0 0 0.243061915040016;0 0 0.242292672395706;0 0 0.241523429751396;0 0 0.240754187107086;0 0 0.239984929561615;0 0 0.239215686917305;0 0 0.19137254357338;0 0 0.143529415130615;0 0 0.0956862717866898;0 0 0.0478431358933449;0 0 0;0.250980406999588 0 0;0.286648005247116 0 0;0.322315603494644 0 0;0.357983201742172 0 0;0.3936507999897 0 0;0.429318398237228 0 0;0.464985996484756 0 0;0.500653624534607 0 0;0.536321222782135 0 0;0.571988821029663 0 0;0.607656419277191 0 0;0.643324017524719 0 0;0.678991615772247 0 0;0.714659214019775 0 0;0.750326812267303 0 0;0.785994410514832 0 0;0.82166200876236 0 0;0.857329607009888 0 0;0.892997205257416 0 0;0.928664803504944 0 0;0.964332401752472 0 0;1 0 0;0.999380826950073 0.042930856347084 0.00990712083876133;0.998761594295502 0.0858617126941681 0.0198142416775227;0.998142421245575 0.128792569041252 0.029721362516284;0.997523248195648 0.171723425388336 0.0396284833550453;0.996904015541077 0.21465428173542 0.0495356060564518;0.99628484249115 0.257585138082504 0.059442725032568;0.995665609836578 0.30051600933075 0.0693498477339745;0.995046436786652 0.343446850776672 0.0792569667100906;0.994427263736725 0.386377722024918 0.0891640856862068;0.993808031082153 0.42930856347084 0.0990712121129036;0.993188858032227 0.472239434719086 0.10897833108902;0.9925696849823 0.515170276165009 0.118885450065136;0.991950452327728 0.558101117610931 0.128792569041252;0.991331279277802 0.601032018661499 0.138699695467949;0.99071204662323 0.643962860107422 0.148606806993484;0.990092873573303 0.686893701553345 0.158513933420181;0.989473700523376 0.729824542999268 0.168421059846878;0.988854467868805 0.772755444049835 0.178328171372414;0.988235294818878 0.815686285495758 0.18823529779911;0.988388061523438 0.818079948425293 0.195670992136002;0.988540887832642 0.820473670959473 0.203106701374054;0.988693654537201 0.822867333889008 0.210542395710945;0.98884642124176 0.825260996818542 0.217978104948997;0.988999247550964 0.827654719352722 0.225413799285889;0.989152014255524 0.830048382282257 0.232849508523941;0.989304840564728 0.832442104816437 0.240285202860832;0.989457607269287 0.834835767745972 0.247720912098885;0.989610373973846 0.837229430675507 0.255156606435776;0.989763200283051 0.839623153209686 0.262592315673828;0.98991596698761 0.842016816139221 0.27002802491188;0.990068733692169 0.844410479068756 0.27746370434761;0.990221560001373 0.846804201602936 0.284899413585663;0.990374326705933 0.849197864532471 0.292335122823715;0.990527093410492 0.851591527462006 0.299770832061768;0.990679919719696 0.853985249996185 0.307206511497498;0.990832686424255 0.85637891292572 0.31464222073555;0.990985512733459 0.8587726354599 0.322077929973602;0.991138279438019 0.861166298389435 0.329513639211655;0.991291046142578 0.86355996131897 0.336949318647385;0.991443872451782 0.865953683853149 0.344385027885437;0.991596639156342 0.868347346782684 0.351820737123489;0.991749405860901 0.870741009712219 0.359256446361542;0.991902232170105 0.873134732246399 0.366692125797272;0.992054998874664 0.875528395175934 0.374127835035324;0.992207765579224 0.877922058105469 0.381563544273376;0.992360591888428 0.880315780639648 0.388999253511429;0.992513358592987 0.882709443569183 0.396434932947159;0.992666184902191 0.885103166103363 0.403870642185211;0.99281895160675 0.887496829032898 0.411306351423264;0.99297171831131 0.889890491962433 0.418742060661316;0.993124544620514 0.892284214496613 0.426177740097046;0.993277311325073 0.894677877426147 0.433613449335098;0.993430078029633 0.897071540355682 0.441049158573151;0.993582904338837 0.899465262889862 0.448484867811203;0.993735671043396 0.901858925819397 0.455920547246933;0.993888437747955 0.904252588748932 0.463356256484985;0.994041264057159 0.906646311283112 0.470791965723038;0.994194030761719 0.909039974212646 0.47822767496109;0.994346857070923 0.911433696746826 0.48566335439682;0.994499623775482 0.913827359676361 0.493099063634872;0.994652390480042 0.916221022605896 0.500534772872925;0.994805216789246 0.918614745140076 0.507970452308655;0.994957983493805 0.921008408069611 0.51540619134903;0.995110750198364 0.923402070999146 0.52284187078476;0.995263576507568 0.925795793533325 0.53027755022049;0.995416343212128 0.92818945646286 0.537713289260864;0.995569109916687 0.930583119392395 0.545148968696594;0.995721936225891 0.932976841926575 0.552584707736969;0.99587470293045 0.93537050485611 0.560020387172699;0.996027529239655 0.937764227390289 0.567456066608429;0.996180295944214 0.940157890319824 0.574891805648804;0.996333062648773 0.942551553249359 0.582327485084534;0.996485888957977 0.944945275783539 0.589763164520264;0.996638655662537 0.947338938713074 0.597198903560638;0.996791422367096 0.949732601642609 0.604634582996368;0.9969442486763 0.952126324176788 0.612070322036743;0.997097015380859 0.954519987106323 0.619506001472473;0.997249782085419 0.956913650035858 0.626941680908203;0.997402608394623 0.959307372570038 0.634377419948578;0.997555375099182 0.961701035499573 0.641813099384308;0.997708201408386 0.964094758033752 0.649248778820038;0.997860968112946 0.966488420963287 0.656684517860413;0.998013734817505 0.968882083892822 0.664120197296143;0.998166561126709 0.971275806427002 0.671555936336517;0.998319327831268 0.973669469356537 0.678991615772247;0.998472094535828 0.976063132286072 0.686427295207977;0.998624920845032 0.978456854820251 0.693863034248352;0.998777687549591 0.980850517749786 0.701298713684082;0.99893045425415 0.983244180679321 0.708734393119812;0.999083280563354 0.985637903213501 0.716170132160187;0.999236047267914 0.988031566143036 0.723605811595917;0.999388873577118 0.990425288677216 0.731041550636292;0.999541640281677 0.99281895160675 0.738477230072021;0.999694406986237 0.995212614536285 0.745912909507751;0.999847233295441 0.997606337070465 0.753348648548126;1 1 0.760784327983856];

imagesc(Brain_roi_weight05n)
caxis([-0.015 0.018])
colormap (map)
for j = 1:length(indexb) % draw lines dividing network
    if j<2
        line_number =indexb(j)-0.1;
	else
        line_number = indexb(j)-0.01; 
    end
    line([0,355],[line_number,line_number],'Color','white','Linewidth',2);%画网络的分割线
    line([line_number,line_number],[0,355],'Color','white','Linewidth',2);
end
ax = gca;
name_order2 = network_name_new;
tick_order2 = index;
set(ax, 'XTick',tick_order2 , 'XTickLabel', name_order2);
xtickangle(90)
set(ax, 'YTick',tick_order2, 'YTickLabel', name_order2);
set(ax,'FontSize',10,'FontName','Times New Roman','FontWeight','bold');
box off
axis square;
width=5000;
height=5000;
saveas(gcf,'D:\projects\ABCD\EF\PLS\results\ef_task23\combat\2fold\comp2\roi_weight.png')





