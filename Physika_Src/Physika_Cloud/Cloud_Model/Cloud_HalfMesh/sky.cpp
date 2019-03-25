#include "sky.h"

Sky::Sky()
{
	cout << "creat a sky object" << endl;
}

void Sky::Initial()
{
	img_sky = NULL;
	img_grey_sky = NULL;
}

void Sky::CreateSkyPossion(Image image, Pixel pixel)//���ò�ֵ��������ճ�����������ջҶȵ�ֵ���ڸ߶ȳ����㹫ʽ
{
	int*  isSkyPixelList = new int[image.GetImg_width()*image.GetImg_height()];
	for (int i = 0; i<image.GetImg_height(); i++)
		for (int j = 0; j<image.GetImg_width(); j++)
		{
			isSkyPixelList[i*image.GetImg_width() + j] = pixel.GetPixelTypeList()[i*image.GetImg_width() + j];
		}


	int neighbor_radius = 5;
	for (int i = 0; i<image.GetImg_height(); i++)
		for (int j = 0; j<image.GetImg_width(); j++)
		{
			bool isFind = false;
			for (int k = -neighbor_radius; k<neighbor_radius && !isFind; k++)
				for (int n = -neighbor_radius; n<neighbor_radius && !isFind; n++)
				{
					int idx = j + n;
					int idy = i + k;
					if (idx>0 && idx<image.GetImg_width() &&idy>0 && idy<image.GetImg_height())
					{
						if (pixel.GetPixelTypeList()[idy*image.GetImg_width() + idx] != 0)
						{
							isFind = true;
						}
					}

				}

			if (isFind)
			{
				isSkyPixelList[i*image.GetImg_width() + j] = 0;
			}

			else
			{
				isSkyPixelList[i*image.GetImg_width() + j] = 1;
			}


		}


	img_sky = new Color4f[image.GetImg_width()*image.GetImg_height()];
	Color4f*  tempCor = new Color4f[image.GetImg_width()*image.GetImg_height()];

	Color4f avgSky(0, 0, 0);		//average sky color		
	int skycount = 0;
	int cloudPixelCount = 0;
	for (int i = 0; i<image.GetImg_width()*image.GetImg_height(); i++)
	{
		if (isSkyPixelList[i] == 1)
		{
			img_sky[i] = image.GetImg_mat_cor()[i];

			//avgSky += image.GetImg_mat_cor()[i];

			avgSky.setRedChannel(avgSky.redChannel() + image.GetImg_mat_cor()[i].redChannel());
			avgSky.setGreenChannel(avgSky.greenChannel() + image.GetImg_mat_cor()[i].greenChannel());
			avgSky.setBlueChannel(avgSky.blueChannel() + image.GetImg_mat_cor()[i].blueChannel());

			skycount++;
		}


	}
	//avgSky = avgSky*(1.0 / skycount);

	avgSky.setRedChannel(avgSky.redChannel()*(1.0 / skycount));
	avgSky.setGreenChannel(avgSky.greenChannel()*(1.0 / skycount));
	avgSky.setBlueChannel(avgSky.blueChannel()*(1.0 / skycount));

	avgSky = Color4f(0.02, 0.02, 0.02);

	//set colors of cloud  pixel to avgsky
	cloudPixelCount = 0;
	for (int i = 0; i<image.GetImg_height(); i++)
		for (int j = 0; j<image.GetImg_width(); j++)
		{
			if (isSkyPixelList[i*image.GetImg_width() + j] == 1)
			{
				img_sky[i*image.GetImg_width() + j] = avgSky;
				cloudPixelCount++;

			}
		}

	int Max_Iter = 5000;

	float Error = 999999999;

	neighbor_radius = 3;

	//�ر�ʱ�䣡��
	for (int iter = 0; iter<Max_Iter&&Error>1.0e-20; iter++)
	{
		for (int i = 0; i<image.GetImg_width()*image.GetImg_height(); i++)
		{
			tempCor[i] = Color4f(0, 0, 0);
		}

		for (int i = 0; i<image.GetImg_height(); i++)
			for (int j = 0; j<image.GetImg_width(); j++)
			{
				if (isSkyPixelList[i*image.GetImg_width() + j] == 0 && i >= 0 && i<image.GetImg_height()&&j >= 0 && j<image.GetImg_width())
				{
					int neighorcount = 0;
					float weightSum = 0.0;
					for (int k = -neighbor_radius; k<neighbor_radius; k++)
						for (int n = -neighbor_radius; n<neighbor_radius; n++)
						{
							if (n*n + k*k != 0)
							{
								int idx = j + n;
								int idy = i + k;
								if (idx>0 && idx<image.GetImg_width()&&idy>0 && idy<image.GetImg_height())
								{
									if (isSkyPixelList[idy*image.GetImg_width() + idx] == 1)
									{
										//tempCor[i*image.GetImg_width() + j] += image.GetImg_mat_cor()[idy*image.GetImg_width() + idx];

										tempCor[i*image.GetImg_width() + j].setRedChannel(tempCor[i*image.GetImg_width() + j].redChannel() + image.GetImg_mat_cor()[idy*image.GetImg_width() + idx].redChannel());
										tempCor[i*image.GetImg_width() + j].setGreenChannel(tempCor[i*image.GetImg_width() + j].greenChannel() + image.GetImg_mat_cor()[idy*image.GetImg_width() + idx].greenChannel());
										tempCor[i*image.GetImg_width() + j].setBlueChannel(tempCor[i*image.GetImg_width() + j].blueChannel() + image.GetImg_mat_cor()[idy*image.GetImg_width() + idx].blueChannel());
									}
										
									else
									{
										//tempCor[i*image.GetImg_width() + j] += img_sky[idy*image.GetImg_width() + idx];

										tempCor[i*image.GetImg_width() + j].setRedChannel(tempCor[i*image.GetImg_width() + j].redChannel() + img_sky[idy*image.GetImg_width() + idx].redChannel());
										tempCor[i*image.GetImg_width() + j].setGreenChannel(tempCor[i*image.GetImg_width() + j].greenChannel() + img_sky[idy*image.GetImg_width() + idx].greenChannel());
										tempCor[i*image.GetImg_width() + j].setBlueChannel(tempCor[i*image.GetImg_width() + j].blueChannel() + img_sky[idy*image.GetImg_width() + idx].blueChannel());
									}
										

									neighorcount++;

								}
							}
						}

					//Color4f temp = tempCor[i*image.GetImg_width() + j] * (1.0 / neighorcount);

					Color4f temp;
					temp.setRedChannel(tempCor[i*image.GetImg_width() + j].redChannel()*(1.0 / neighorcount));
					temp.setGreenChannel(tempCor[i*image.GetImg_width() + j].greenChannel()*(1.0 / neighorcount));
					temp.setBlueChannel(tempCor[i*image.GetImg_width() + j].blueChannel()*(1.0 / neighorcount));

					//tempCor[i*image.GetImg_width() + j] = tempCor[i*image.GetImg_width() + j] * (1.0 / neighorcount);

					tempCor[i*image.GetImg_width() + j].setRedChannel(tempCor[i*image.GetImg_width() + j].redChannel()*(1.0 / neighorcount));
					tempCor[i*image.GetImg_width() + j].setGreenChannel(tempCor[i*image.GetImg_width() + j].greenChannel()*(1.0 / neighorcount));
					tempCor[i*image.GetImg_width() + j].setBlueChannel(tempCor[i*image.GetImg_width() + j].blueChannel()*(1.0 / neighorcount));
				}

			}

		float curError = 0;
		for (int i = 0; i<image.GetImg_width()*image.GetImg_height(); i++)
		{
			if (isSkyPixelList[i] == 0)
			{
				//Color4f difCor = tempCor[i] - img_sky[i];

				Color4f difCor;
				difCor.setRedChannel(tempCor[i].redChannel() - img_sky[i].redChannel());
				difCor.setGreenChannel(tempCor[i].greenChannel() - img_sky[i].greenChannel());
				difCor.setBlueChannel(tempCor[i].blueChannel() - img_sky[i].blueChannel());

				curError += difCor.redChannel()*difCor.redChannel() + difCor.greenChannel()*difCor.greenChannel() + difCor.blueChannel()*difCor.blueChannel();



			}
		}

		curError = sqrt(curError / cloudPixelCount);


		cout<<"IterCount:" <<iter<<"curError: "<<curError<<endl;
		for (int i = 0; i<image.GetImg_width()*image.GetImg_height(); i++)
		{
			if (isSkyPixelList[i] == 1)
			{
				img_sky[i] = image.GetImg_mat_cor()[i];

			}
			else
				img_sky[i] = tempCor[i];

		}
		Error = curError;

	}


	//convert the colorful  sky to grey image
	img_grey_sky = new float[image.GetImg_width()*image.GetImg_height()];

	for (int i = 0; i<image.GetImg_height(); i++)
		for (int j = 0; j<image.GetImg_width(); j++)
		{
			float R = img_sky[i*image.GetImg_width() + j].redChannel();
			float G = img_sky[i*image.GetImg_width() + j].greenChannel();
			float B = img_sky[i*image.GetImg_width() + j].blueChannel();
			img_grey_sky[i*image.GetImg_width() + j] = 0.2989 * R + 0.5870 * G + 0.1140 * B;//0.2989 * R + 0.5870 * G + 0.1140 * B ;

		}


	float max_sky_dif_cloud = -99999;

	int nSingularPixel = 0;
	for (int i = 0; i<image.GetImg_height(); i++)
		for (int j = 0; j<image.GetImg_width(); j++)
		{
			if (pixel.GetPixelTypeList()[i*image.GetImg_width() + j] == 1)
			{
				if (img_grey_sky[i*image.GetImg_width() + j] - image.GetImg_mat()[i*image.GetImg_width() + j]>0)
				{
					nSingularPixel++;

				}

			}
		}
	cout << "The number of pixels whose sky intensities being greater than image intensities:   " << nSingularPixel << endl;

	delete[] isSkyPixelList;
}

float* Sky::GetImg_grey_sky()
{
	return this->img_grey_sky;
}
