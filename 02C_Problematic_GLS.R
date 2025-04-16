
removed_trips <- c(
  '1620006001', '1646001001',  '6101002001', 
  '13218001003', '17275002001', '19922003001', '19959003001', '2214006001',
  '2293005003', '23311004001', '23319003001', '23315003001',
  '23334002001', '23351003001', '23352004001', '23355004001', 'BD763004001',
  'BD763005001', 'BD807004001', 'BM258001001', 'BR666001001', 'BX322001001',
  'J606002001',
  'S941001001', 'T036002001', 'V396015002001'
)

moved_trips <- c(
  '1643001001', '19426003001', '19426003002',  '19955003001',
  '21166002001', '23297004001', 
  '23299003001', 
  '23348004001',  'N647003001', 'S977002001',
  'T024002001', 'T049001001', 'T049002001', 'T064001001', 'T064002001', 
  'V396018002001'
)

cut_trips <- c(
  '1605002003', '1619005001', '1647005001', '1711006001', '23291004001',
  '23307004001', '23315004001', '23335004001','BD887004001', 'J549002001', 
  'J576002001', 'J594002001',
  'N799002001', 'N810001002', 'N813002001', 'N817002001'
)

incomplete_trips <- c(
  '1598002001', '1600002001', '1603003001', '1605001001', '1611001001',
  '1613001001', '1613002001', '1613003001', '1617002001', '1619002001',
  '1620001001', '1622002001', '1622003001', '1629001001', '1630003001',
  '1638002004', '1640001001', '1647001001', '1648002001', '1652001001',
  '1681002001', '1703001001', '1711002001', '2280003001', '8364001001',
  '8599002001', '19088001002', '19205002001', '19211001001', '23290004001',
  '23293001001', '23300001001', '23348001001', '23362003001', 'BD768004001',
  'BD800003001', 'BD808001001', 'BD811001001', 'BD812004001', 'BR50328001001',
  'BR50347001002', 'BX317001001', 'BX339001001', 'J484002001', 'J550002001',
  'J572002001', 'J579002001', 'J580001002', 'J585001002', 'J588001001', 
  'J589001001', 'N208004001', 'N788001001', 'N792002001', 'N795003001',
  'N811003001', 'P961002001', 'P965002001', 'Q666002001', 'T045003001',
  'T053002001', 'V396028001001', 'V396034002001', 'X236001001', 'X242001001',
  'X247001001', 'X250001001', 'X255001001', 'X465001001'
)

removed <- data.frame(geo_trip = removed_trips,
                    Problem = 'Removed')
moved <- data.frame(geo_trip = moved_trips,
                  Problem = 'Moved')
cut <- data.frame(geo_trip = cut_trips,
                  Problem = 'Cut')
incomplete <- data.frame(geo_trip = incomplete_trips,
                         Problem = 'Incomplete')
bad_trips <- bind_rows (removed, moved, cut, incomplete)

trips <- read_csv('D:/Dropbox/Diego Vicente/TESIS/01_GLS_analysis/Results/SGAT/Raw_migratory_trips_tol.csv')
cut_dates <- read_csv('D:/Dropbox/Diego Vicente/TESIS/01_GLS_analysis/Raw Data/cut_dates_resuscitated.csv', col_names = TRUE)
# gls_removed <- read_csv('D:/Dropbox/Diego Vicente/TESIS/01_GLS_analysis/Raw Data/removed_resuscitated.csv', col_names = TRUE)


ids <- data.frame(geo_trip = unique(trips$geo_trip),
                  i = c(1:length(unique(trips$geo_trip))))

bad_trips <- bad_trips %>%
  left_join(cut_dates)%>%
  left_join(ids)

write_csv(bad_trips, 'D:/Dropbox/Diego Vicente/TESIS/01_GLS_analysis/Raw Data/bad_trips.csv')
                  