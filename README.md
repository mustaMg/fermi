# Sonrasında Uzun Süreli Işıma Evresi Gösteren Gama Işını Patlamalarının Bulunması
Bu repoda TUBİTAK STAR projesi kapsamında yaptığım EE arayan kodlar bulunmakta.

[English](https://github.com/mustafagumustas/Search-of-Extended-Emissions-on-GRB-s/blob/main/README_en.md)

## Gama Işın Patlamaları (GRB) nedir?
Samanyolu Galaksisi dışından gelen çok yüksek enerjilere sahip, kozmik olaylardır. GRB’ler öncül ışınımdan gelen akı sürelerine göre uzun veya kısa süreli olmak üzere iki kategoriye ayrılırlar, bu süreye T90 denir. T90, dedektöre gelen fotonların %90’ının ne kadar sürede toplandığıdır. Uzun GRB’ler T90>2 sn ve yumuşak spektrumlu, kısa olanlar ise T90<2 sn ve sert spektrumlu olurlar. GRB’lerin spektrumlarının sertlik ve yumuşaklığı yüksek enerji kanalındaki toplam akısının düşük enerji kanalındaki toplam akısına oranıyla belirlenir, yüksek orana sahip olanlar (>1) sert spektrumludurlar.

  İncelenen bazı GRB’ler yukarıda tanımlanan kategorilere net olarak uymamışlardır. Öncül ışınımı incelendiğinde kısa olarak tanımlanması gereken GRB’ler, bu öncül ışınımı takip eden uzun süreli ışımaları (extended emission: EE) göstermesi nedeniyle uzun patlama olarak kategorilen- dirilmiştir. Her iki kategorinin de özelliğini gösterebilen bu hibrit kategorinin yani EE’li GRB’lerin incelenmesi uzun ve kısa GRB’lerin daha iyi anlaşılmasını sağlayacaktır.

## EE Araması
  EE’leri incelemek adına yapılan bu çalışma için 10 adet GRB seçildi. Swift tarafından oluşturulan kayıtlı verilerden “Işık Eğrisi” (Light Curve, LC) oluşturularak morfolojik kriterlere göre seçim yapıldı. LC’leri oluşturmak için Swift Portal platformu üzerinden GRB’lerin ‘log’, ‘bat’ ve ‘auxil’ dosyaları indirildi. Bat verileri içerisindeki event dosyaları ‘HEASOFT’ içerisinde tanımlı Swift komutları ile 5 çeşit LC dosyası üretildi. Bunlar:
- Tek enerji kanallı (15-150 eV), 1 saniye çözünürlüklü, zemin ışıması çıkartılmış,
- Tek enerji kanallı, 64 milisaniye çözünürlüklü, zemin ışıması çıkartılmış,
- 4 enerji kanallı (15-25, 25-50, 50-100, 100-150 eV), 1 saniye çözünürlüklü, zemin ışıması
çıkartılmış,
- 4 enerji kanallı, 64 milisaniye çözünürlüklü, zemin ışıması çıkartılmış LC dosyalarıdır.
- 4 enerji kanallı, 1 saniye çözünürlüklü, filtre uygulanmamış (unweighted) LC dosyalarıdır.

## Morfolojik Kriterlere Göre Verilerin Filtrelenmesi
T90 süreleri 5sn'den uzun olan bazı patlamalar, öncü ışınım sonrasında uzun süreli ışımaya (EE) sahip olabilir. Bunu araştırmak için morfolojik kriterler şu şekilde tanımlanmıştır: GRB’nin maksimum foton sayısına ulaştığı an ilk 5 saniye içerisinde bulunmalı. Bu maksimum değer 11000’den (foton sayısı/ saniye) küçük ise %40, büyük ise %30’una denk gelen ilk andan 5. saniyeye kadar olan sürenin %60’ında; foton sayımları bu değerin altında olması gerekir.

## Kod
Yukarıda anlatılan morfolojik kriterler python ile kodlaranrak fits dosyalarının okunup içerisinde EE olan bölge varsa işaretleyerek pdf biçiminde kaydedildi.
